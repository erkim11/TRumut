/*
catalog.cpp:
  - k-mer catalog loading
  - 31-mer encoding and scanning
  - output rows
*/

#include <zlib.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>

#include <algorithm>
#include <string>
#include <unordered_set>
#include <mutex>
#include <vector>

#include "catalog.h"
#include "third_party/kom-alloc.h"
#include "third_party/kseq.h"

KSTREAM_INIT(gzFile, gzread, 0x10000)

// encode a 31-mer into a 2-bit key; -1 if it contains a non-ACGT base.
int encode_31mer(const char *seq31, uint64_t *out_key)
{
    uint64_t x = 0;
    for (int i = 0; i < 31; ++i) {
        uint8_t c = kom_nt4_table[(uint8_t)seq31[i]];
        if (c >= 4) return -1;
        x = (x << 2) | (uint64_t)c;
    }
    *out_key = x;
    return 0;
}

int split_tab(char *s, char **f, int max_f)
{
    int n = 0;
    while (*s && n < max_f) {
        f[n++] = s;
        char *t = strchr(s, '\t');
        if (!t) break;
        *t = 0;
        s = t + 1;
    }
    return n;
}

const char *detect_header(int with_qname, int with_coords, int with_seq)
{
    static char buf[256];
    snprintf(buf, sizeof buf,
             "chr\tregion_start\tregion_end\tloci\t%s"
             "kmerA_id\tkmerB_id\t%s"
             "n_reads\tallele_lengths\tn_motif_kmers\t%sHP\tPS\n",
             with_qname ? "query_names\t" : "",
             with_coords ? "kmerA_ref_coords\tkmerB_ref_coords\t" : "",
             with_seq ? "allele_seqs\t" : "");
    return buf;
}

typedef std::unordered_set<uint64_t> MotifSet;

// per-pair, per-output accumulator.
struct PairAcc {
    uint32_t n = 0;
    std::string lens, motifs, hp, ps, names, seqs;
};

void pair_free_extras(pair_t *p)
{
    delete (MotifSet *)p->motif_kmers;
    p->motif_kmers = 0;
    for (int w = 0; w < 2; ++w) { delete (PairAcc *)p->acc[w]; p->acc[w] = 0; }
}

static inline void list_append_cstr(std::string &list, uint32_t n, const char *s)
{
    if (n) list += ',';
    list += s;
}

static const MotifSet *pair_motif_kmers(pair_t *pair)
{
    if (!pair->motif_kmers)
        pair->motif_kmers = new MotifSet(find_loci_motif_kmers(pair->loci ? pair->loci : ""));
    return (const MotifSet *)pair->motif_kmers;
}

void pair_precompute_motifs(pair_t *pairs, int32_t n_pairs)
{
    for (int32_t i = 0; i < n_pairs; ++i) pair_motif_kmers(&pairs[i]);
}

enum { ACC_STRIPES = 4096 };
static std::mutex acc_stripes[ACC_STRIPES];

// record one read's hit
void record_hit(pair_t *pair, int32_t pair_index, int which,
                const char *name,
                const token_t *tok, int32_t t1_ref, int32_t t2_ref,
                char kmerA_match, char kmerB_match,
                int32_t posA, int32_t posB,
                const char *seq, int with_seq,
                const char *hp, const char *ps)
{
    int32_t left  = posA < posB ? posA : posB;
    int32_t right = posA < posB ? posB : posA;
    int32_t dist = right - (left + 31);
    if (dist < 0) return;
    const char *between_ptr = seq + left + 31;

    const MotifSet *motifs = pair_motif_kmers(pair);
    unsigned long long n_motif = motifs->empty() ? 0 : count_motif_kmers(between_ptr, dist, *motifs);

    std::lock_guard<std::mutex> lk(acc_stripes[(uint32_t)pair_index % ACC_STRIPES]);

    if (!pair->acc[which]) pair->acc[which] = new PairAcc();
    PairAcc &a = *(PairAcc *)pair->acc[which];
    char buf[32];

    snprintf(buf, sizeof buf, "%d", dist);
    list_append_cstr(a.lens, a.n, buf);

    if (motifs->empty()) list_append_cstr(a.motifs, a.n, ".");
    else {
        snprintf(buf, sizeof buf, "%llu", n_motif);
        list_append_cstr(a.motifs, a.n, buf);
    }

    const int has_hp = hp && *hp && strcmp(hp, "NA") != 0 && strcmp(hp, ".") != 0;
    list_append_cstr(a.hp, a.n, has_hp ? hp : ".");
    const int has_ps = has_hp && ps && *ps && strcmp(ps, "NA") != 0 && strcmp(ps, ".") != 0;
    list_append_cstr(a.ps, a.n, has_ps ? ps : ".");

    if (name) list_append_cstr(a.names, a.n, *name ? name : ".");   // --save-query-names

    if (with_seq) {
        char kmerA_strand = tok[t1_ref].strand;
        char kmerB_strand = tok[t2_ref].strand;
        int antisense = (kmerA_strand != kmerA_match) && (kmerB_strand != kmerB_match);
        if (a.n) a.seqs += ',';
        if (dist == 0) {
            a.seqs += '.';
        } else if (antisense) {
            size_t at = a.seqs.size();
            a.seqs.append(between_ptr, (size_t)dist);
            kom_revcomp((uint64_t)dist, &a.seqs[at]);
        } else {
            a.seqs.append(between_ptr, (size_t)dist);
        }
    }
    a.n++;
}

int flush_pair_rows(FILE *out, pair_t *pairs, int32_t n_pairs, const token_t *tok,
                    int which, int with_qname, int with_coords, int with_seq)
{
    std::vector<int32_t> order;
    order.reserve(1024);
    for (int32_t i = 0; i < n_pairs; ++i)
        if (pairs[i].acc[which] && ((PairAcc *)pairs[i].acc[which])->n) order.push_back(i);

    std::sort(order.begin(), order.end(), [&](int32_t x, int32_t y) {
        const pair_t &p = pairs[x], &q = pairs[y];
        if (p.chr_rank != q.chr_rank) return p.chr_rank < q.chr_rank;
        int c = strcmp(p.chr ? p.chr : "", q.chr ? q.chr : "");
        if (c) return c < 0;
        c = strcmp(tok[p.tok1].id, tok[q.tok1].id);
        if (c) return c < 0;
        return strcmp(tok[p.tok2].id, tok[q.tok2].id) < 0;
    });

    for (int32_t i : order) {
        pair_t &p = pairs[i];
        const PairAcc &a = *(const PairAcc *)p.acc[which];
        if (fprintf(out, "%s\t%lld\t%lld\t%s\t",
                    p.chr ? p.chr : "", (long long)p.region_start, (long long)p.region_end,
                    p.loci ? p.loci : "") < 0) return -1;
        if (with_qname && fprintf(out, "%s\t", a.names.c_str()) < 0) return -1;
        if (fprintf(out, "%s\t%s\t", tok[p.tok1].id, tok[p.tok2].id) < 0) return -1;
        if (with_coords && fprintf(out, "%ld\t%ld\t", (long)p.pos1, (long)p.pos2) < 0) return -1;
        if (fprintf(out, "%u\t%s\t%s\t",
                    (unsigned)a.n, a.lens.c_str(), a.motifs.c_str()) < 0) return -1;
        if (with_seq && fprintf(out, "%s\t", a.seqs.c_str()) < 0) return -1;
        if (fprintf(out, "%s\t%s\n", a.hp.c_str(), a.ps.c_str()) < 0) return -1;
    }
    return fflush(out) == 0 ? 0 : -1;
}

// roll a 31-mer key along the sequence and look each one up in the dictionary, appending hit positions to pv.
void scan_read_for_tokens(const char *read_seq, int32_t read_len,
                          key2i_t *key2i,
                          posvec_t *pv, int32_t **hit_tok, int32_t *n_hit_tok, int32_t *m_hit_tok)
{
    int i, l;
    const int k = 31;
    uint64_t x = 0, mask = (1ULL<<(k*2)) - 1;

    for (i = l = 0, x = 0; i < read_len; ++i) {
        int c = kom_nt4_table[(uint8_t)read_seq[i]];
        if (c < 4) {
            x = (x << 2 | (uint64_t)c) & mask;
            if (++l >= k) {
                int32_t pos = i - k + 1;
                khint_t it = key2i_get(key2i, x);
                if (it != kh_end(key2i)) {
                    int32_t tok = kh_val(key2i, it);
                    if (hit_tok && pv[tok].n == 0) {
                        if (*n_hit_tok == *m_hit_tok) kom_grow(int32_t, *hit_tok, *n_hit_tok, *m_hit_tok);
                        (*hit_tok)[(*n_hit_tok)++] = tok;
                    }
                    if (pv[tok].n == pv[tok].m) kom_grow(int32_t, pv[tok].a, pv[tok].n, pv[tok].m);
                    pv[tok].a[pv[tok].n++] = pos;
                }
            }
        } else {
            l = 0; x = 0;
        }
    }
}

static int32_t get_or_add_token(token_t **tok, int32_t *n_tok, int32_t *m_tok,
                                str2i_t *token_index,
                                const char *id, char strand, const char *seq31)
{
    char keybuf[4096];
    snprintf(keybuf, sizeof(keybuf), "%s\t%c", id, strand);
    khint_t it = str2i_get(token_index, keybuf);
    if (it != kh_end(token_index)) return kh_val(token_index, it);

    int32_t idx = *n_tok;
    if (*n_tok == *m_tok) kom_grow(token_t, *tok, *n_tok, *m_tok);
    (*tok)[idx].id = kom_strdup(id);
    (*tok)[idx].strand = strand;
    for (int i = 0; i < 31; ++i)
        (*tok)[idx].seq[i] = (char)toupper((unsigned char)seq31[i]);
    (*tok)[idx].seq[31] = 0;
    (*tok)[idx].key = 0;

    int absent;
    khint_t it2 = str2i_put(token_index, keybuf, &absent);
    kh_key(token_index, it2) = kom_strdup(keybuf);
    kh_val(token_index, it2) = idx;
    ++(*n_tok);
    return idx;
}

void add_strand_twins(token_t **tok, int32_t *n_tok, int32_t *m_tok,
                      str2i_t *token_index)
{
    int32_t n_orig = *n_tok;
    for (int32_t i = 0; i < n_orig; ++i) {
        char other = (*tok)[i].strand == '+' ? '-' : '+';
        char *id = (*tok)[i].id;
        char seqbuf[32];
        memcpy(seqbuf, (*tok)[i].seq, 32);
        get_or_add_token(tok, n_tok, m_tok, token_index, id, other, seqbuf);
    }
}

static void record_skip(skipped_row_t **skip, int32_t *n_skip, int32_t *m_skip,
                        const char *idA, const char *idB, const char *reason)
{
    if (!skip) return;
    if (*n_skip == *m_skip) kom_grow(skipped_row_t, *skip, *n_skip, *m_skip);
    (*skip)[*n_skip].idA = kom_strdup(idA && *idA ? idA : "NA");
    (*skip)[*n_skip].idB = kom_strdup(idB && *idB ? idB : "NA");
    (*skip)[*n_skip].reason = reason;
    ++(*n_skip);
}

// one "start-end:motif" locus token.
static char *make_locus_token(const char *start, const char *end, const char *motif)
{
    size_t n = strlen(start) + strlen(end) + strlen(motif) + 3;
    char *t = kom_malloc(char, n);
    snprintf(t, n, "%s-%s:%s", start, end, motif);
    return t;
}

static void pair_add_locus(pair_t *p, const char *token)
{
    size_t tl = strlen(token);
    if (p->loci) {
        for (const char *s = p->loci; (s = strstr(s, token)); s += 1) {
            int at_start = (s == p->loci || s[-1] == ',');
            int at_end   = (s[tl] == 0 || s[tl] == ',');
            if (at_start && at_end) return;
        }
        size_t ol = strlen(p->loci);
        char *merged = kom_malloc(char, ol + 1 + tl + 1);
        memcpy(merged, p->loci, ol);
        merged[ol] = ',';
        memcpy(merged + ol + 1, token, tl + 1);
        free(p->loci);
        p->loci = merged;
    } else {
        p->loci = kom_strdup(token);
    }
}

int chr_name_eq(const char *a, const char *b)
{
    if (strncasecmp(a, "chr", 3) == 0) a += 3;
    if (strncasecmp(b, "chr", 3) == 0) b += 3;
    return strcasecmp(a, b) == 0;
}

int read_pairs_file(const char *pairs_path, const char *chr_filter,
                    pair_t **pairs, int32_t *n_pairs, int32_t *m_pairs,
                    token_t **tok, int32_t *n_tok, int32_t *m_tok,
                    str2i_t *token_index, // key = "id\tstrand" -> tok idx
                    skipped_row_t **skip, int32_t *n_skip, int32_t *m_skip)
{
    gzFile fp = gzopen(pairs_path, "rb");
    if (!fp) return -1;

    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};

    // (tok1, tok2) -> pairs index, for this file's rows.
    key2i_t *pair_index = key2i_init();

    enum { MAX_FIELDS = 1024 };
    char *fields[MAX_FIELDS];
    int64_t n_bad_seq = 0;

    static const char *const names[] = {
        "chr", "tr_start", "tr_end", "tr_motif",
        "kmerA_ref_coords", "kmerA_id", "kmerA_seq", "kmerA_strand",
        "kmerB_ref_coords", "kmerB_id", "kmerB_seq", "kmerB_strand" };
    enum { C_CHR, C_START, C_END, C_MOTIF, C_A_POS, C_A_ID, C_A_SEQ, C_A_STRAND,
           C_B_POS, C_B_ID, C_B_SEQ, C_B_STRAND, N_NAMES };
    int col[N_NAMES];
    int n_cols = 0;
    {
        int got = 0;
        while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0)
            if (str.l > 0) { got = 1; break; }
        if (!got) {
            fprintf(stderr, "[ERROR] k-mer catalog %s is empty\n", pairs_path);
            key2i_destroy(pair_index); free(str.s); ks_destroy(ks); gzclose(fp);
            return -2;
        }
        char *hdr = kom_strdup(str.s);
        int nh = split_tab(hdr, fields, MAX_FIELDS);
        if (nh > 0) { size_t l = strlen(fields[nh-1]); if (l && fields[nh-1][l-1] == '\r') fields[nh-1][l-1] = 0; }
        for (int c = 0; c < N_NAMES; ++c) {
            col[c] = -1;
            for (int i = 0; i < nh; ++i)
                if (strcmp(fields[i], names[c]) == 0) { col[c] = i; break; }
            if (col[c] < 0) {
                fprintf(stderr, "[ERROR] %s: catalog header lacks column '%s' (rebuild with trcall prep-catalog)\n",
                        pairs_path, names[c]);
                free(hdr); key2i_destroy(pair_index); free(str.s); ks_destroy(ks); gzclose(fp);
                return -2;
            }
            if (col[c] + 1 > n_cols) n_cols = col[c] + 1;
        }
        free(hdr);
    }

    while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
        if (str.l == 0) continue;
        char *line = kom_strdup(str.s);

        int nf = split_tab(line, fields, MAX_FIELDS);

        if (chr_filter && nf > col[C_CHR] && !chr_name_eq(fields[col[C_CHR]], chr_filter)) {
            free(line); continue;
        }

        if (nf < n_cols) {
            record_skip(skip, n_skip, m_skip,
                        nf > col[C_A_ID] ? fields[col[C_A_ID]] : 0,
                        nf > col[C_B_ID] ? fields[col[C_B_ID]] : 0,
                        "fewer columns than the header");
            free(line); continue;
        }

        int64_t pos1 = strtoll(fields[col[C_A_POS]], 0, 10);
        const char *id1 = fields[col[C_A_ID]];
        char *seq1 = fields[col[C_A_SEQ]];
        char strand1 = fields[col[C_A_STRAND]][0];

        int64_t pos2 = strtoll(fields[col[C_B_POS]], 0, 10);
        const char *id2 = fields[col[C_B_ID]];
        char *seq2 = fields[col[C_B_SEQ]];
        char strand2 = fields[col[C_B_STRAND]][0];

        if (!((strand1 == '+' || strand1 == '-') && (strand2 == '+' || strand2 == '-'))) {
            record_skip(skip, n_skip, m_skip, id1, id2, "invalid strand column");
            free(line); continue;
        }

        int badA = strlen(seq1) != 31, badB = strlen(seq2) != 31;
        if (badA || badB) {
            ++n_bad_seq;
            record_skip(skip, n_skip, m_skip, id1, id2,
                        badA && badB ? "kmerA_seq and kmerB_seq not 31-mers"
                        : badA       ? "kmerA_seq not a 31-mer"
                                     : "kmerB_seq not a 31-mer");
            free(line); continue;
        }

        if (strand1 == '-') kom_revcomp(31, seq1);
        if (strand2 == '-') kom_revcomp(31, seq2);

        int32_t t1 = get_or_add_token(tok, n_tok, m_tok, token_index, id1, strand1, seq1);
        int32_t t2 = get_or_add_token(tok, n_tok, m_tok, token_index, id2, strand2, seq2);

        int64_t str_start = strtoll(fields[col[C_START]], 0, 10);
        int64_t str_end   = strtoll(fields[col[C_END]], 0, 10);
        char *locus = make_locus_token(fields[col[C_START]], fields[col[C_END]], fields[col[C_MOTIF]]);

        uint64_t pkey = ((uint64_t)(uint32_t)t1 << 32) | (uint32_t)t2;
        int absent;
        khint_t it = key2i_put(pair_index, pkey, &absent);
        if (absent) {
            if (*n_pairs == *m_pairs) kom_grow(pair_t, *pairs, *n_pairs, *m_pairs);
            pair_t *p = &(*pairs)[*n_pairs];
            p->tok1 = t1;
            p->tok2 = t2;
            p->pos1 = pos1;
            p->pos2 = pos2;
            p->chr  = kom_strdup(fields[col[C_CHR]]);
            p->chr_rank = chr_rank_cstr(p->chr);
            p->region_start = str_start;
            p->region_end   = str_end;
            p->loci = 0;
            p->motif_kmers = 0;
            p->acc[0] = p->acc[1] = 0;
            pair_add_locus(p, locus);
            kh_val(pair_index, it) = *n_pairs;
            ++(*n_pairs);
        } else {
            pair_t *p = &(*pairs)[kh_val(pair_index, it)];
            if (str_start < p->region_start) p->region_start = str_start;
            if (str_end   > p->region_end)   p->region_end   = str_end;
            pair_add_locus(p, locus);
        }
        free(locus);

        free(line);
    }

    if (n_bad_seq)
        fprintf(stderr, "[WARNING] %s: skipped %lld row(s) with a non-31-mer kmerA_seq/kmerB_seq\n", pairs_path, (long long)n_bad_seq);

    key2i_destroy(pair_index);
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return 0;
}

// encode each token's 31-mer (reverse-complemented for '-' tokens) and index it as key -> token.
int build_dictionary(token_t *tok, int32_t n_tok, key2i_t *key2i)
{
    for (int32_t ti = 0; ti < n_tok; ++ti) {
        char seq31[32];
        memcpy(seq31, tok[ti].seq, 32);

        if (tok[ti].strand == '-')
            kom_revcomp(31, seq31);

        uint64_t key;
        if (encode_31mer(seq31, &key) != 0) {
            fprintf(stderr, "[ERROR] kmer %s contains non-ACGT after strand handling\n", tok[ti].id);
            return -1;
        }
        tok[ti].key = key;

        int absent;
        khint_t it_put = key2i_put(key2i, key, &absent);
        kh_val(key2i, it_put) = ti;
    }
    return 0;
}
