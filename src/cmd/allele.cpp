/*
allele.cpp -- `trumut allele`: calls tandem repeat alleles, their sequences, and their lengths
*/
#include <zlib.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <stddef.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <condition_variable>
#include <deque>
#include <mutex>
#include <thread>
#include <vector>

#include "lib/catalog.h"
#include "lib/input.h"
#include "third_party/ketopt.h"
#include "third_party/khashl.h"
#include "third_party/kom-alloc.h"
#include "third_party/kseq.h"

KSEQ_INIT2(static __attribute__((unused)), gzFile, gzread)

typedef struct {
    pair_t   *pairs;   int32_t n_pairs;
    token_t  *tok;     int32_t n_tok;
    key2i_t  *key2i;
    str2i_t  *token_index;
    skipped_row_t *skip; int32_t n_skip;   /* catalog rows that never loaded */
} dict_t;

typedef struct {
    const char *input, *ref, *target_chr;
    const char *catalog;
    const char *out_path, *error_path, *invert_path, *debug_path;
    int mapq_min, n_threads, fmt;
    int n_jobs;
    int shard;
    int save_qnames;
    int save_seq;
    int save_coords;
    int ploidy;
    int reads_mode;
    int is_reads;
    int inflate_threads;
    seqsrc_t *pre_src;
} opt_t;

static int path_exists(const char *p)
{
    struct stat st;
    return stat(p, &st) == 0;
}

/* distinct names in the catalog's chr column, in first-appearance order */
static int catalog_chr_list(const char *path, char ***out)
{
    *out = 0;
    gzFile fp = gzopen(path, "rb");
    if (!fp) return -1;
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};

    char **names = 0;
    int n = 0, m = 0;
    int chr_col = -1;
    while (chr_col < 0 && ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
        if (str.l == 0) continue;
        int i = 0;
        for (char *p = str.s, *e = str.s + str.l; p < e; ++i) {
            char *t = (char *)memchr(p, '\t', (size_t)(e - p));
            size_t l = t ? (size_t)(t - p) : (size_t)(e - p);
            if (l && p[l-1] == '\r') --l;
            if (l == 3 && strncmp(p, "chr", 3) == 0) { chr_col = i; break; }
            if (!t) break;
            p = t + 1;
        }
        if (chr_col < 0) {
            fprintf(stderr, "[ERROR] k-mer catalog %s: header line has no 'chr' column\n", path);
            free(str.s); ks_destroy(ks); gzclose(fp);
            return -1;
        }
    }
    while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
        if (str.l == 0) continue;
        char *p = str.s, *e = str.s + str.l;
        for (int i = 0; i < chr_col && p; ++i) {
            char *t = (char *)memchr(p, '\t', (size_t)(e - p));
            p = t ? t + 1 : 0;
        }
        if (!p) continue;
        char *tab = (char *)memchr(p, '\t', (size_t)(e - p));
        size_t len = tab ? (size_t)(tab - p) : (size_t)(e - p);
        if (len == 0) continue;
        memmove(str.s, p, len);

        if (n > 0 && strlen(names[n-1]) == len && strncmp(names[n-1], str.s, len) == 0)
            continue;
        int seen = 0;
        for (int i = 0; i < n && !seen; ++i)
            seen = strlen(names[i]) == len && strncmp(names[i], str.s, len) == 0;
        if (seen) continue;

        if (n == m) { m = m ? m * 2 : 32; names = (char **)realloc(names, (size_t)m * sizeof *names); }
        names[n] = (char *)malloc(len + 1);
        memcpy(names[n], str.s, len);
        names[n][len] = 0;
        ++n;
    }
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    *out = names;
    return n;
}

static int find_index(const char *path, const char *ext, const char *base_ext,
                      char *out, size_t n)
{
    snprintf(out, n, "%s%s", path, ext);
    if (path_exists(out)) return 0;

    size_t pl = strlen(path), bl = strlen(base_ext);
    if (pl > bl && strcasecmp(path + pl - bl, base_ext) == 0) {
        snprintf(out, n, "%.*s%s", (int)(pl - bl), path, ext);
        if (path_exists(out)) return 0;
    }
    return -1;
}

static int require_index(const opt_t *o)
{
    if (o->fmt != SR_FMT_BAM && o->fmt != SR_FMT_CRAM) return 0;

    const char *ext  = o->fmt == SR_FMT_BAM ? ".bai"  : ".crai";
    const char *base = o->fmt == SR_FMT_BAM ? ".bam"  : ".cram";
    const char *kind = sr_fmt_name(o->fmt);

    if (strcmp(o->input, "-") == 0) {
        if (o->fmt == SR_FMT_BAM) return 0;
        fprintf(stderr, "[ERROR] %s input needs an index and cannot be read from stdin; pass the %s file path\n", kind, kind);
        return -1;
    }

    char idx[PATH_MAX];
    if (find_index(o->input, ext, base, idx, sizeof idx) != 0) {
        fprintf(stderr, "[ERROR] no %s index for %s (run: samtools index %s); looked for %s%s and %s -> %s\n",
                ext, o->input, o->input, o->input, ext, base, ext);
        return -1;
    }

    struct stat sd, si;
    if (stat(o->input, &sd) == 0 && stat(idx, &si) == 0 && si.st_mtime < sd.st_mtime)
        fprintf(stderr, "[INFO] %s is older than %s; re-run samtools index if the data changed\n",
                idx, o->input);

    return 0;
}

#ifdef __APPLE__
extern "C" int sysctlbyname(const char *, void *, size_t *, void *, size_t);
#endif

static int parse_int_opt(const char *flag, const char *s, int *out)
{
    char *end = 0;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (!s[0] || !end || *end || errno || v < INT_MIN || v > INT_MAX) {
        fprintf(stderr, "[ERROR] %s must be an integer (got '%s')\n", flag, s);
        return -1;
    }
    *out = (int)v;
    return 0;
}

static int cpu_count(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int)n : 1;
}

static long long usable_ram_bytes(void)
{
#ifdef __APPLE__
    /* Darwin has no cheap MemAvailable, so use total RAM, discounted below. */
    uint64_t mem = 0;
    size_t len = sizeof mem;
    if (sysctlbyname("hw.memsize", &mem, &len, 0, 0) == 0 && mem > 0)
        return (long long)mem;
    return 0;
#else
    FILE *f = fopen("/proc/meminfo", "r");
    if (f) {
        char line[256];
        long long kb = 0;
        while (fgets(line, sizeof line, f)) {
            if (sscanf(line, "MemAvailable: %lld kB", &kb) == 1) { fclose(f); return kb * 1024; }
        }
        fclose(f);
    }
    long pages = sysconf(_SC_AVPHYS_PAGES), psz = sysconf(_SC_PAGESIZE);
    if (pages > 0 && psz > 0) return (long long)pages * psz;
    pages = sysconf(_SC_PHYS_PAGES);
    if (pages > 0 && psz > 0) return (long long)pages * psz;
    return 0;
#endif
}
#define CATALOG_MEM_FACTOR 18
#define WORKER_BASE_BYTES  (256LL << 20)

static long long worker_mem_estimate(const char *catalog, int n_shards)
{
    struct stat st;
    if (stat(catalog, &st) != 0 || st.st_size <= 0) return 0;
    long long share = (long long)st.st_size * 2 / (n_shards > 0 ? n_shards : 1);
    if (share > st.st_size) share = st.st_size;
    return share * CATALOG_MEM_FACTOR + WORKER_BASE_BYTES;
}

#define DEFAULT_JOBS 3   /* -j default everywhere: small and fixed, never the whole machine */

static int auto_jobs(long long per_worker, int n_shards, char *why, size_t why_n)
{
    int cores = cpu_count();
    int j = DEFAULT_JOBS;
    const char *limit = "the default of 3";

    if (j > n_shards)      { j = n_shards;      limit = "chromosomes to run"; }

    long long ram = usable_ram_bytes();
    if (ram > 0 && per_worker > 0) {
        long long fit = (long long)((double)ram * 0.8) / per_worker;
        if (fit < 1) fit = 1;
        if (fit < j) { j = (int)fit; limit = "memory"; }
    }
    if (j < 1) j = 1;

    snprintf(why, why_n, "%d cores, %.1f GB usable, ~%.0f MB per worker; limited by %s",
             cores, ram / 1073741824.0, per_worker / 1048576.0, limit);
    return j;
}

/* columns of the -m (--many-kmer-hits) file. */
#define DETECT_ERROR_HEADER \
    "query_name\tchr\tregion_start\tregion_end\tkmerA_id\tkmerB_id\tkmerA_hit_count\tkmerB_hit_count\n"

/* columns of the -d file. */
#define DETECT_DEBUG_HEADER "kmerA_id\tkmerB_id\tdropped\treason\n"

static inline uint64_t pack_tt(int32_t t1, int32_t t2) {
    return ((uint64_t)(uint32_t)t1 << 32) | (uint32_t)t2;
}

static int open_debug(const opt_t *o, const dict_t *d, gzFile *dbg)
{
    *dbg = 0;
    if (!o->debug_path || !o->debug_path[0]) return 0;
    *dbg = strcmp(o->debug_path, "-") == 0
         ? gzdopen(dup(STDOUT_FILENO), "wb")
         : gzopen(o->debug_path, "wb");
    if (!*dbg) { fprintf(stderr, "[ERROR] cannot open debug file %s\n", o->debug_path); return -1; }
    if (!o->shard) gzputs(*dbg, DETECT_DEBUG_HEADER);

    str2i_t *seen = str2i_init();
    char keybuf[4096];
    int absent;
    for (int32_t i = 0; i < d->n_pairs; ++i) {
        snprintf(keybuf, sizeof(keybuf), "%s\t%s",
                 d->tok[d->pairs[i].tok1].id, d->tok[d->pairs[i].tok2].id);
        khint_t it = str2i_put(seen, keybuf, &absent);
        if (absent) kh_key(seen, it) = kom_strdup(keybuf);
    }
    for (int32_t i = 0; i < d->n_skip; ++i) {
        snprintf(keybuf, sizeof(keybuf), "%s\t%s", d->skip[i].idA, d->skip[i].idB);
        if (str2i_get(seen, keybuf) != kh_end(seen)) continue;
        snprintf(keybuf, sizeof(keybuf), "%s\t%s\t%s",
                 d->skip[i].idA, d->skip[i].idB, d->skip[i].reason);
        khint_t it = str2i_put(seen, keybuf, &absent);
        if (!absent) continue;
        kh_key(seen, it) = kom_strdup(keybuf);
        gzprintf(*dbg, "%s\t%s\tyes\tnot loaded: %s\n",
                 d->skip[i].idA, d->skip[i].idB, d->skip[i].reason);
    }
    for (khint_t it = 0; it != kh_end(seen); ++it)
        if (kh_exist(seen, it)) free((char *)kh_key(seen, it));
    str2i_destroy(seen);
    return 0;
}

static const char *header_for(const opt_t *o)
{
    return detect_header(o->save_qnames, o->save_coords, o->save_seq);
}

static int wants_gz(const char *path)
{
    if (!path || strcmp(path, "-") == 0) return 0;
    size_t n = strlen(path);
    return n > 3 && strcmp(path + n - 3, ".gz") == 0;
}

static int write_header_file(const char *path, const char *header)
{
    if (wants_gz(path)) {
        gzFile g = gzopen(path, "wb");
        if (!g) return -1;
        int rc = (header && gzputs(g, header) < 0) ? -1 : 0;
        if (gzclose(g) != Z_OK) rc = -1;
        return rc;
    }
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    if (header) fputs(header, f);
    return fclose(f) == 0 ? 0 : -1;
}

static FILE *open_rows_output(const char *path)
{
    if (!path) return stdout;
    if (wants_gz(path)) {
        char cmd[PATH_MAX + 64];
        snprintf(cmd, sizeof cmd, "gzip -c >> '%s'", path);
        fflush(NULL);
        return popen(cmd, "w");
    }
    return fopen(path, "ab");
}

static int close_rows_output(FILE *f, const char *path)
{
    if (!f) return 0;
    if (f == stdout) return fflush(f) == 0 ? 0 : -1;
    if (path && wants_gz(path)) {
        int st = pclose(f);
        if (st == -1 || !WIFEXITED(st) || WEXITSTATUS(st) != 0) {
            fprintf(stderr, "[ERROR] the output compressor failed\n");
            return -1;
        }
        return 0;
    }
    return fclose(f) == 0 ? 0 : -1;
}

static int finish_outputs(const opt_t *o, pair_t *pairs, int32_t n_pairs, const token_t *tok,
                          FILE *out, FILE *invert_out)
{
    const char *out_file = (o->out_path && strcmp(o->out_path, "-") != 0) ? o->out_path : 0;
    int rc = 0;
    if (flush_pair_rows(out, pairs, n_pairs, tok, 0, o->save_qnames, o->save_coords, o->save_seq) != 0) {
        fprintf(stderr, "[ERROR] writing the output failed\n"); rc = -1;
    }
    if (close_rows_output(out, out_file) != 0) rc = -1;
    if (invert_out) {
        if (flush_pair_rows(invert_out, pairs, n_pairs, tok, 1, o->save_qnames, o->save_coords, o->save_seq) != 0) {
            fprintf(stderr, "[ERROR] writing the invert file failed\n"); rc = -1;
        }
        if (close_rows_output(invert_out, o->invert_path) != 0) rc = -1;
    }
    return rc;
}

static int open_outputs(const opt_t *o, FILE **out, FILE **invert_out)
{
    const char *out_file = (o->out_path && strcmp(o->out_path, "-") != 0) ? o->out_path : 0;
    if (out_file) {
        if (write_header_file(out_file, o->shard ? 0 : header_for(o)) != 0) {
            fprintf(stderr, "[ERROR] cannot open output %s\n", o->out_path); return -1;
        }
    } else if (!o->shard) {
        fputs(header_for(o), stdout);
    }
    *out = open_rows_output(out_file);
    if (!*out) { fprintf(stderr, "[ERROR] cannot open output %s\n", out_file ? out_file : "-"); return -1; }

    *invert_out = 0;
    if (o->invert_path && o->invert_path[0]) {
        if (write_header_file(o->invert_path, o->shard ? 0 : header_for(o)) != 0) {
            fprintf(stderr, "[ERROR] cannot open invert file %s\n", o->invert_path); return -1;
        }
        *invert_out = open_rows_output(o->invert_path);
        if (!*invert_out) { fprintf(stderr, "[ERROR] cannot open invert file %s\n", o->invert_path); return -1; }
    }
    return 0;
}

enum { PS_EMIT = 1, PS_INVERT = 2, PS_OVERLAP = 4, PS_MULTI = 8 };

typedef struct {
    const dict_t *d;
    const opt_t *o;
    int32_t *partner;            /* partner[t]: the same k-mer id on the opposite strand */
    token_pairs_t *token2pairs;  /* token -> the pairs it appears in */
    int want_invert;             /* -w given: record strand/order mismatches too */
} scan_shared_t;

typedef struct {
    posvec_t *pv;                /* pv[tok]: that token's positions in the current read */
    int32_t *hit_tok;            /* the tokens hit, so only those need resetting */
    int32_t n_hit_tok, m_hit_tok;
    u64set_t *seen_pairs;        /* pairs already handled in the current read */
    uint8_t *pair_state;         /* -d: fates per pair, else NULL */
    uint8_t *tok_seen;           /* -d: tokens ever hit, else NULL */
    FILE *err_out;               /* -e rows, or NULL */
} scan_local_t;

static int scan_shared_init(scan_shared_t *s, const dict_t *d, const opt_t *o)
{
    const pair_t *pairs = d->pairs; int32_t n_pairs = d->n_pairs;
    const token_t *tok = d->tok;   int32_t n_tok = d->n_tok;

    memset(s, 0, sizeof *s);
    s->d = d;
    s->o = o;
    s->want_invert = o->invert_path && o->invert_path[0];

    // partner[t]: the same k-mer id on the opposite strand, for revcomp reads
    s->partner = (int32_t *)kom_calloc(int32_t, (size_t)n_tok);
    {
        char keybuf[4096];
        for (int32_t t = 0; t < n_tok; ++t) {
            char other = tok[t].strand == '+' ? '-' : '+';
            snprintf(keybuf, sizeof(keybuf), "%s\t%c", tok[t].id, other);
            khint_t it = str2i_get(d->token_index, keybuf);
            s->partner[t] = (it != kh_end(d->token_index)) ? (int32_t)kh_val(d->token_index, it) : -1;
        }
    }

    // token -> the pairs it appears in
    token_pairs_t *token2pairs = kom_calloc(token_pairs_t, n_tok);
    for (int32_t pi = 0; pi < n_pairs; ++pi) {
        int32_t t1 = pairs[pi].tok1;
        int32_t t2 = pairs[pi].tok2;

        if (token2pairs[t1].n == token2pairs[t1].m)
            kom_grow(int32_t, token2pairs[t1].pairs, token2pairs[t1].n, token2pairs[t1].m);
        token2pairs[t1].pairs[token2pairs[t1].n++] = pi;

        if (t2 != t1) {   // skip a self-pair
            if (token2pairs[t2].n == token2pairs[t2].m)
                kom_grow(int32_t, token2pairs[t2].pairs, token2pairs[t2].n, token2pairs[t2].m);
            token2pairs[t2].pairs[token2pairs[t2].n++] = pi;
        }
    }
    s->token2pairs = token2pairs;
    return 0;
}

static void scan_shared_free(scan_shared_t *s)
{
    if (s->token2pairs)
        for (int32_t i = 0; i < s->d->n_tok; ++i) free(s->token2pairs[i].pairs);
    free(s->token2pairs);
    free(s->partner);
    memset(s, 0, sizeof *s);
}

static int scan_local_init(scan_local_t *c, const scan_shared_t *s, FILE *err_out, int want_debug)
{
    memset(c, 0, sizeof *c);
    c->pv = kom_calloc(posvec_t, s->d->n_tok);
    c->seen_pairs = u64set_init();
    if (want_debug) {
        c->pair_state = kom_calloc(uint8_t, (size_t)s->d->n_pairs);
        c->tok_seen   = kom_calloc(uint8_t, (size_t)s->d->n_tok);
    }
    c->err_out = err_out;
    return 0;
}

static void scan_local_free(scan_local_t *c, const scan_shared_t *s)
{
    if (c->pv) for (int32_t i = 0; i < s->d->n_tok; ++i) free(c->pv[i].a);
    free(c->pv);
    free(c->hit_tok);
    if (c->seen_pairs) u64set_destroy(c->seen_pairs);
    free(c->pair_state);
    free(c->tok_seen);
    memset(c, 0, sizeof *c);
}

/* scan one read and record every catalog pair it spans. */
static void process_read(const scan_shared_t *s, scan_local_t *c,
                         const char *qname, const char *read_seq, int32_t L,
                         const char *hp, const char *ps)
{
    const dict_t *d = s->d;
    const opt_t *o = s->o;
    pair_t *pairs = d->pairs;
    const token_t *tok = d->tok;
    const int32_t *partner = s->partner;
    const token_pairs_t *token2pairs = s->token2pairs;
    posvec_t *pv = c->pv;

    c->n_hit_tok = 0;
    scan_read_for_tokens(read_seq, L, d->key2i, pv, &c->hit_tok, &c->n_hit_tok, &c->m_hit_tok);

    u64set_destroy(c->seen_pairs);
    c->seen_pairs = u64set_init();
    u64set_t *seen_pairs = c->seen_pairs;

    for (int32_t k = 0; k < c->n_hit_tok; ++k) {
        int32_t tok_idx = c->hit_tok[k];
        if (pv[tok_idx].n == 0) continue;

        int32_t src[2] = { tok_idx, partner[tok_idx] };
        for (int sidx = 0; sidx < 2; ++sidx) {
        if (src[sidx] < 0 || (sidx == 1 && src[1] == src[0])) continue;

        for (int32_t j = 0; j < token2pairs[src[sidx]].n; ++j) {
            int32_t pi = token2pairs[src[sidx]].pairs[j];

            int absent;
            u64set_put(seen_pairs, (uint64_t)pi, &absent);
            if (!absent) continue;

            int32_t t1 = pairs[pi].tok1;
            int32_t t2 = pairs[pi].tok2;
            int32_t t1_alt = partner[t1] >= 0 ? partner[t1] : t1;
            int32_t t2_alt = partner[t2] >= 0 ? partner[t2] : t2;
            int32_t na1 = pv[t1].n, na2 = pv[t1_alt].n;
            int32_t nb1 = pv[t2].n, nb2 = pv[t2_alt].n;
            int32_t nA = na1 + na2;
            int32_t nB = nb1 + nb2;
            if (nA == 0 || nB == 0) continue;

            if (nA != 1 || nB != 1) {
                if (c->pair_state) c->pair_state[pi] |= PS_MULTI;
                if (c->err_out) {
                    const char *chr_s = pairs[pi].chr ? pairs[pi].chr : "";
                    fprintf(c->err_out, "%s\t%s\t%lld\t%lld\t%s\t%s\t%d\t%d\n",
                            qname, chr_s,
                            (long long)pairs[pi].region_start,
                            (long long)pairs[pi].region_end,
                            tok[t1].id, tok[t2].id, nA, nB);
                }
                continue;
            }

            // exactly one hit each: check the distance, then emit once.
            int32_t posA = (na1 > 0) ? pv[t1].a[0] : pv[t1_alt].a[0];
            int32_t posB = (nb1 > 0) ? pv[t2].a[0] : pv[t2_alt].a[0];
            int32_t left = posA < posB ? posA : posB;
            int32_t right = posA < posB ? posB : posA;
            if (right - (left + 31) < 0) { if (c->pair_state) c->pair_state[pi] |= PS_OVERLAP; continue; }
            // a row goes to -w when the two k-mers disagree about strand, or agree but sit in the wrong order for it.
            // *_match is relative to the strand the pair records: '+' matched the recorded 31-mer, '-' matched its reverse complement.
            char kmerA_strand = tok[t1].strand, kmerB_strand = tok[t2].strand;
            char kmerA_match = (na1 > 0) ? '+' : '-';
            char kmerB_match = (nb1 > 0) ? '+' : '-';
            int xor_mismatch = (kmerA_strand != kmerA_match) != (kmerB_strand != kmerB_match);
            int both_match = (kmerA_strand == kmerA_match && kmerB_strand == kmerB_match);
            int both_mismatch = (kmerA_strand != kmerA_match && kmerB_strand != kmerB_match);
            int to_invert = xor_mismatch
                || (both_match && posA > posB)
                || (both_mismatch && posB > posA);
            if (to_invert) {
                if (c->pair_state) c->pair_state[pi] |= PS_INVERT;
                if (s->want_invert)
                    record_hit(&pairs[pi], pi, 1, o->save_qnames ? qname : 0, tok, t1, t2, kmerA_match, kmerB_match, posA, posB, read_seq, o->save_seq, hp, ps);
            } else {
                if (c->pair_state) c->pair_state[pi] |= PS_EMIT;
                record_hit(&pairs[pi], pi, 0, o->save_qnames ? qname : 0, tok, t1, t2, kmerA_match, kmerB_match, posA, posB, read_seq, o->save_seq, hp, ps);
            }
        }
        }
    }

    for (int32_t k = 0; k < c->n_hit_tok; ++k) {
        int32_t t = c->hit_tok[k];
        if (c->tok_seen) c->tok_seen[t] = 1;
        pv[t].n = 0;
    }
}

static void write_debug_classification(gzFile dbg_out, const scan_shared_t *s,
                                       const uint8_t *pair_state, const uint8_t *tok_seen)
{
    const dict_t *d = s->d;
    const pair_t *pairs = d->pairs;
    const token_t *tok = d->tok;
    const int32_t *partner = s->partner;

    u64set_t *dbg_seen = u64set_init();
    for (int32_t pi = 0; pi < d->n_pairs; ++pi) {
        int32_t t1 = pairs[pi].tok1, t2 = pairs[pi].tok2;
        int absent;
        u64set_put(dbg_seen, pack_tt(t1, t2), &absent);
        if (!absent) continue;
        uint8_t st = pair_state[pi];
        const char *reason = 0;
        if      (st & PS_EMIT)    reason = 0;
        else if (st & PS_INVERT)  reason = "strand/order mismatch";
        else if (st & PS_OVERLAP) reason = "kmers overlap";
        else if (st & PS_MULTI)   reason = "kmer(s) hit more than once in a read";
        else {
            int seenA = tok_seen[t1] || (partner[t1] >= 0 && tok_seen[partner[t1]]);
            int seenB = tok_seen[t2] || (partner[t2] >= 0 && tok_seen[partner[t2]]);
            if      (seenA && seenB) reason = "kmers never found in the same read";
            else if (seenA)          reason = "kmerB not found in any read";
            else if (seenB)          reason = "kmerA not found in any read";
            else                     reason = "neither kmer found in any read";
        }
        gzprintf(dbg_out, "%s\t%s\t%s\t%s\n", tok[t1].id, tok[t2].id,
                 reason ? "yes" : "no", reason ? reason : "NA");
    }
    u64set_destroy(dbg_seen);
}

static int run_alignment(const dict_t *d, const opt_t *o)
{
    sr_t *rd = sr_open(o->input, o->ref, o->target_chr, o->n_threads);
    if (!rd) { fprintf(stderr, "[ERROR] %s\n", sr_open_error()); return 2; }

    if (o->target_chr && sr_has_ref(rd, o->target_chr) == 0) {
        fprintf(stderr, "[ERROR] chromosome %s not in the %s header\n", o->target_chr, sr_format(rd));
        sr_close(rd);
        return 2;
    }

    FILE *out = 0, *invert_out = 0;
    if (open_outputs(o, &out, &invert_out) != 0) { sr_close(rd); return 2; }

    FILE *err_out = 0;
    if (o->error_path && o->error_path[0]) {
        err_out = fopen(o->error_path, "wb");
        if (!err_out) { fprintf(stderr, "[ERROR] cannot open --many-kmer-hits file %s\n", o->error_path); sr_close(rd); return 2; }
        if (!o->shard) fputs(DETECT_ERROR_HEADER, err_out);
    }

    gzFile dbg_out = 0;
    if (open_debug(o, d, &dbg_out) != 0) { sr_close(rd); return 2; }

    scan_shared_t shared;
    scan_local_t local;
    scan_shared_init(&shared, d, o);
    scan_local_init(&local, &shared, err_out, dbg_out != 0);
    pair_precompute_motifs(d->pairs, d->n_pairs);

    // reads already seen, by qname hash.
    u64set_t *seen = u64set_init();

    sr_rec_t rec;
    int sr_ret;

    while ((sr_ret = sr_read(rd, &rec)) > 0) {
        // primary alignments only.
        if (rec.flag & (SR_FSECONDARY | SR_FSUPPLEMENTARY)) continue;
        if ((int)rec.mapq < o->mapq_min) continue;

        uint64_t h = kh_hash_str(rec.qname);
        int absent;
        u64set_put(seen, h, &absent);
        if (!absent) continue;

        process_read(&shared, &local, rec.qname, rec.seq, rec.l_seq, rec.hp, rec.ps);
    }

    int rc = 0;
    if (sr_ret < 0) {
        fprintf(stderr, "[ERROR] %s\n", sr_error(rd));
        rc = 2;
    }
    if (sr_close(rd) != 0) rc = 2;

    if (dbg_out) write_debug_classification(dbg_out, &shared, local.pair_state, local.tok_seen);

    u64set_destroy(seen);
    if (finish_outputs(o, d->pairs, d->n_pairs, d->tok, out, invert_out) != 0) rc = 2;
    if (err_out) fclose(err_out);
    if (dbg_out) gzclose(dbg_out);

    scan_local_free(&local, &shared);
    scan_shared_free(&shared);
    return rc;
}

typedef struct {
    uint32_t name_off, seq_off;   /* into read_batch_t::buf */
    uint32_t hp_off, ps_off;      /* 0 when the read carries none */
    int32_t  len;
} read_item_t;

struct read_batch_t {
    std::vector<char> buf;
    std::vector<read_item_t> items;
    int64_t bases = 0;

    uint32_t put(const char *s, size_t n, int upper) {
        uint32_t off = (uint32_t)buf.size();
        buf.insert(buf.end(), s, s + n);
        if (upper)
            for (size_t i = off; i < buf.size(); ++i) buf[i] = (char)toupper((unsigned char)buf[i]);
        buf.push_back(0);
        return off;
    }
    void add(const char *name, const char *seq, int32_t len, const char *hp, const char *ps, int upper) {
        read_item_t it;
        it.name_off = put(name, strlen(name), 0);
        it.seq_off  = put(seq, (size_t)len, upper);
        it.hp_off   = hp ? put(hp, strlen(hp), 0) : 0;
        it.ps_off   = ps ? put(ps, strlen(ps), 0) : 0;
        it.len = len;
        items.push_back(it);
        bases += len;
    }
};

struct batch_queue_t {
    std::mutex mu;
    std::condition_variable cv_put, cv_get;
    std::deque<read_batch_t *> q;
    size_t cap = 1;
    bool closed = false;

    void put(read_batch_t *b) {
        std::unique_lock<std::mutex> lk(mu);
        cv_put.wait(lk, [&] { return q.size() < cap; });
        q.push_back(b);
        cv_get.notify_one();
    }
    read_batch_t *get() {
        std::unique_lock<std::mutex> lk(mu);
        cv_get.wait(lk, [&] { return !q.empty() || closed; });
        if (q.empty()) return 0;
        read_batch_t *b = q.front();
        q.pop_front();
        cv_put.notify_one();
        return b;
    }
    void close() {
        std::lock_guard<std::mutex> lk(mu);
        closed = true;
        cv_get.notify_all();
    }
};

enum { READS_BATCH_N = 256, READS_BATCH_BASES = 4 << 20 };

static int run_reads(const dict_t *d, const opt_t *o)
{
    FILE *out = 0, *invert_out = 0;
    if (open_outputs(o, &out, &invert_out) != 0) return 2;

    FILE *err_out = 0;
    if (o->error_path && o->error_path[0]) {
        err_out = fopen(o->error_path, "wb");
        if (!err_out) { fprintf(stderr, "[ERROR] cannot open --many-kmer-hits file %s\n", o->error_path); return 2; }
        if (!o->shard) fputs(DETECT_ERROR_HEADER, err_out);
    }

    gzFile dbg_out = 0;
    if (open_debug(o, d, &dbg_out) != 0) return 2;

    int n_workers = o->n_jobs > 0 ? o->n_jobs : DEFAULT_JOBS;
    if (n_workers < 1) n_workers = 1;

    seqsrc_t *src = 0;
    sr_t *rd = 0;
    seqrec_t *rr = 0;
    int fmt = o->fmt;
    const char *kind_note = "";
    if (fmt == SR_FMT_URL) {
        rd = sr_open(o->input, o->ref, 0, o->n_threads);
        if (!rd) { fprintf(stderr, "[ERROR] %s\n", sr_open_error()); return 2; }
        kind_note = ", via samtools";
    } else {
        src = o->pre_src;
        if (!src) {
            char err[256];
            src = seqsrc_open(o->input, o->inflate_threads, err, sizeof err);
            if (!src) { fprintf(stderr, "[ERROR] %s\n", err); return 2; }
        }

        if (fmt == SR_FMT_UNKNOWN) {
            uint8_t m[16];
            ssize_t n = seqsrc_peek(src, m, sizeof m);
            if (n <= 0) {
                fprintf(stderr, "[ERROR] %s: %s\n", o->input, n < 0 ? seqsrc_error(src) : "empty input");
                seqsrc_close(src);
                return 2;
            }
            fmt = sr_classify_bytes(m, (size_t)n);
        }
        if (sr_fmt_is_alignment(fmt)) {
            rd = sr_open_source(src, strcmp(o->input, "-") == 0 ? "stdin" : o->input);
            if (!rd) { fprintf(stderr, "[ERROR] %s\n", sr_open_error()); seqsrc_close(src); return 2; }
        } else {
            if (o->target_chr) {
                fprintf(stderr, "[ERROR] -t needs aligned input; %s is %s (unaligned reads carry no chromosome)\n",
                        strcmp(o->input, "-") == 0 ? "stdin" : o->input, sr_fmt_name(fmt));
                seqsrc_close(src);
                return 2;
            }
            rr = seqrec_open(src);
        }
        kind_note = seqsrc_kind(src) == SEQSRC_BGZF
            ? (o->inflate_threads > 1 ? ", BGZF inflated on -@ threads" : ", BGZF")
            : seqsrc_kind(src) == SEQSRC_GZIP ? ", gzip (inflated by one thread, which may bound the run)" : "";
    }
    fprintf(stderr, "[INFO] scanning %s (%s%s) with %d thread%s%s\n",
            rd ? "alignments as one stream" : "unaligned reads",
            rd ? sr_format(rd) : sr_fmt_name(fmt), kind_note, n_workers, n_workers == 1 ? "" : "s",
            o->target_chr ? ", keeping only reads on the -t chromosome" : "");

    scan_shared_t shared;
    scan_shared_init(&shared, d, o);
    pair_precompute_motifs(d->pairs, d->n_pairs);

    batch_queue_t queue;
    queue.cap = (size_t)n_workers * 4 + 2;

    std::vector<scan_local_t> locals((size_t)n_workers);
    for (int i = 0; i < n_workers; ++i) scan_local_init(&locals[i], &shared, err_out, dbg_out != 0);

    std::vector<std::thread> workers;
    for (int i = 0; i < n_workers; ++i) {
        workers.emplace_back([&, i]() {
            scan_local_t *c = &locals[(size_t)i];
            for (;;) {
                read_batch_t *b = queue.get();
                if (!b) break;
                const char *base = b->buf.data();
                for (const read_item_t &it : b->items)
                    process_read(&shared, c, base + it.name_off, base + it.seq_off, it.len,
                                 it.hp_off ? base + it.hp_off : 0, it.ps_off ? base + it.ps_off : 0);
                delete b;
            }
        });
    }

    int rc = 0;
    long long n_reads = 0, n_bases = 0;
    read_batch_t *cur = new read_batch_t();
    auto flush = [&]() {
        if (cur->items.empty()) return;
        queue.put(cur);
        cur = new read_batch_t();
    };
    if (rd) {
        sr_rec_t rec;
        int r;
        while ((r = sr_read(rd, &rec)) > 0) {
            if (rec.flag & (SR_FSECONDARY | SR_FSUPPLEMENTARY)) continue;
            if ((int)rec.mapq < o->mapq_min) continue;
            if (o->target_chr && !chr_name_eq(rec.rname, o->target_chr)) continue;
            cur->add(rec.qname, rec.seq, rec.l_seq, rec.hp, rec.ps, 0);
            ++n_reads; n_bases += rec.l_seq;
            if (cur->items.size() >= READS_BATCH_N || cur->bases >= READS_BATCH_BASES) flush();
        }
        if (r < 0) { fprintf(stderr, "[ERROR] %s\n", sr_error(rd)); rc = 2; }
        if (sr_close(rd) != 0) rc = 2;
    } else {
        const char *name, *seq;
        int32_t len;
        int r;
        while ((r = seqrec_next(rr, &name, &seq, &len)) > 0) {
            cur->add(name, seq, len, 0, 0, o->save_seq);
            ++n_reads; n_bases += len;
            if (cur->items.size() >= READS_BATCH_N || cur->bases >= READS_BATCH_BASES) flush();
        }
        if (r < 0) {
            const char *e = seqsrc_error(src);
            if (*e) fprintf(stderr, "[ERROR] %s: %s\n", o->input, e);
            rc = 2;
        }
        seqrec_close(rr);
    }
    flush();
    delete cur;
    queue.close();
    for (auto &t : workers) t.join();
    if (src && seqsrc_close(src) != 0) rc = 2;

    fprintf(stderr, "[INFO] scanned %lld reads (%lld bases)\n", n_reads, n_bases);

    if (dbg_out) {
        uint8_t *pair_state = kom_calloc(uint8_t, (size_t)d->n_pairs);
        uint8_t *tok_seen   = kom_calloc(uint8_t, (size_t)d->n_tok);
        for (int i = 0; i < n_workers; ++i) {
            for (int32_t p = 0; p < d->n_pairs; ++p) pair_state[p] |= locals[(size_t)i].pair_state[p];
            for (int32_t t = 0; t < d->n_tok; ++t)   tok_seen[t]   |= locals[(size_t)i].tok_seen[t];
        }
        write_debug_classification(dbg_out, &shared, pair_state, tok_seen);
        free(pair_state);
        free(tok_seen);
    }

    if (finish_outputs(o, d->pairs, d->n_pairs, d->tok, out, invert_out) != 0) rc = 2;
    if (err_out) fclose(err_out);
    if (dbg_out) gzclose(dbg_out);

    for (int i = 0; i < n_workers; ++i) scan_local_free(&locals[(size_t)i], &shared);
    scan_shared_free(&shared);
    return rc;
}

#define SEQ_MAX_CONTIGS_CAP 8
static int seq_max_contigs = 2;   /* set from --ploidy before the run */
typedef struct {
    int32_t idx;
    uint64_t seq;
    int32_t count;                              /* total hits, all contigs */
    int32_t n_contigs;                          /* distinct contigs hit, capped at seq_max_contigs + 1 */
    int32_t contig_idx[SEQ_MAX_CONTIGS_CAP];
    int32_t location[SEQ_MAX_CONTIGS_CAP];
    char    match_strand[SEQ_MAX_CONTIGS_CAP];  /* strand of the token that matched */
} kmer_hits_t;

typedef struct {
    char   *seq;
    int32_t len;
    int32_t posA, posB;
} anchor_slice_t;

typedef struct {
    kmer_hits_t A, B;
    anchor_slice_t slice[SEQ_MAX_CONTIGS_CAP];
} pair_hit_info_t;

static int note_kmer_hits(kmer_hits_t *h, int32_t contig_idx,
                          const posvec_t *pv, int32_t t, int32_t t_rc)
{
    int32_t n1 = pv[t].n, n2 = (t_rc >= 0) ? pv[t_rc].n : 0;
    int32_t n = n1 + n2;
    if (n == 0) return -1;
    int slot = -1;
    if (h->n_contigs < seq_max_contigs) {
        int32_t pos1 = (n1 > 0) ? pv[t].a[0] : INT32_MAX;
        int32_t pos2 = (n2 > 0) ? pv[t_rc].a[0] : INT32_MAX;
        int s = h->n_contigs;
        if (pos1 <= pos2) { h->location[s] = pv[t].a[0];    h->match_strand[s] = '+'; }
        else              { h->location[s] = pv[t_rc].a[0]; h->match_strand[s] = '-'; }
        h->contig_idx[s] = contig_idx;
        slot = s;
    }
    if (h->n_contigs <= seq_max_contigs) h->n_contigs++;
    h->count += n;
    return slot;
}

static const char *pair_hits_problem(const pair_hit_info_t *info)
{
    const kmer_hits_t *A = &info->A, *B = &info->B;
    if (A->count == 0 && B->count == 0) return "neither kmer found";
    if (A->count == 0) return "kmerA not found";
    if (B->count == 0) return "kmerB not found";
    int multiA = A->n_contigs > seq_max_contigs || A->count != A->n_contigs;
    int multiB = B->n_contigs > seq_max_contigs || B->count != B->n_contigs;
    if (multiA && multiB) return "both kmers found multiple times";
    if (multiA) return "kmerA found multiple times";
    if (multiB) return "kmerB found multiple times";
    if (A->n_contigs != B->n_contigs) return "kmers on different contigs";
    for (int i = 0; i < A->n_contigs; ++i) {
        int found = 0;
        for (int j = 0; j < B->n_contigs; ++j) found |= (B->contig_idx[j] == A->contig_idx[i]);
        if (!found) return "kmers on different contigs";
    }
    return 0;
}

static int pair_hits_partner(const pair_hit_info_t *info, int i)
{
    for (int j = 0; j < info->B.n_contigs && j < seq_max_contigs; ++j)
        if (info->B.contig_idx[j] == info->A.contig_idx[i]) return j;
    return -1;
}

KHASHL_MAP_INIT(KH_LOCAL, contig2idx_t, contig2idx, const char*, int32_t, kh_hash_str, kh_eq_str)

static int32_t add_contig(char ***names, int32_t *n, int32_t *m, const char *name)
{
    if (*n == *m) {
        *m = *n + 1;
        *m += (*m >> 1) + 16;
        *names = (char**)kom_realloc(char*, *names, *m);
    }
    (*names)[*n] = kom_strdup(name ? name : "");
    return (*n)++;
}

static void stash_slice(pair_hit_info_t *info, int s, const char *seq, int32_t seq_len,
                        int32_t posA, int32_t posB)
{
    int32_t left = posA < posB ? posA : posB;
    int32_t right = posA < posB ? posB : posA;
    int32_t end = right + 31;
    if (end > seq_len) end = seq_len;
    int32_t len = end - left;
    if (len < 0) len = 0;
    anchor_slice_t *sl = &info->slice[s];
    free(sl->seq);
    sl->seq = kom_malloc(char, len + 1);
    for (int32_t i = 0; i < len; ++i) sl->seq[i] = (char)toupper((unsigned char)seq[left + i]);
    sl->seq[len] = 0;
    sl->len = len;
    sl->posA = posA - left;
    sl->posB = posB - left;
}

/* scan one contig. */
static void scan_contig(const scan_shared_t *sh, scan_local_t *c,
                        pair_hit_info_t *hits, int32_t contig_idx,
                        const char *seq, int32_t seq_len)
{
    const dict_t *d = sh->d;
    const int32_t *partner = sh->partner;
    const token_pairs_t *token2pairs = sh->token2pairs;
    posvec_t *pv = c->pv;

    u64set_destroy(c->seen_pairs);
    c->seen_pairs = u64set_init();

    for (int32_t k = 0; k < c->n_hit_tok; ++k) {
        int32_t tok_idx = c->hit_tok[k];
        if (pv[tok_idx].n == 0) continue;
        int32_t src[2] = { tok_idx, partner[tok_idx] };
        for (int sidx = 0; sidx < 2; ++sidx) {
            if (src[sidx] < 0 || (sidx == 1 && src[1] == src[0])) continue;
            for (int32_t j = 0; j < token2pairs[src[sidx]].n; ++j) {
                int32_t pi = token2pairs[src[sidx]].pairs[j];
                int absent;
                u64set_put(c->seen_pairs, (uint64_t)pi, &absent);
                if (!absent) continue;

                int32_t t1 = d->pairs[pi].tok1, t2 = d->pairs[pi].tok2;
                pair_hit_info_t *info = &hits[pi];
                int sa = note_kmer_hits(&info->A, contig_idx, pv, t1, partner[t1]);
                int sb = note_kmer_hits(&info->B, contig_idx, pv, t2, partner[t2]);
                if (sa >= 0 && sb >= 0)
                    stash_slice(info, sa, seq, seq_len, info->A.location[sa], info->B.location[sb]);
            }
        }
    }
}

static int run_sequence(const dict_t *d, const opt_t *o)
{
    pair_t *pairs = d->pairs; int32_t n_pairs = d->n_pairs;
    token_t *tok = d->tok;
    const char *target_chr = o->target_chr;

    pair_hit_info_t *hits = kom_calloc(pair_hit_info_t, (size_t)n_pairs);
    for (int32_t pi = 0; pi < n_pairs; ++pi) {
        pair_hit_info_t *info = &hits[pi];
        info->A.idx = pairs[pi].tok1; info->A.seq = tok[pairs[pi].tok1].key;
        info->B.idx = pairs[pi].tok2; info->B.seq = tok[pairs[pi].tok2].key;
        for (int s = 0; s < SEQ_MAX_CONTIGS_CAP; ++s) {
            info->A.contig_idx[s] = info->B.contig_idx[s] = -1;
            info->A.location[s] = info->B.location[s] = -1;
        }
    }
    char **contig_names = 0;
    int32_t n_contigs = 0, m_contigs = 0;

    scan_shared_t shared;
    scan_shared_init(&shared, d, o);
    scan_local_t local;
    scan_local_init(&local, &shared, 0, 0);

    seqsrc_t *src = o->pre_src;
    if (!src) {
        char err[256];
        src = seqsrc_open(o->input, o->inflate_threads, err, sizeof err);
        if (!src) { fprintf(stderr, "[ERROR] %s\n", err); return 2; }
    }
    seqrec_t *rr = seqrec_open(src);
    int rc = 0;
    int contig_count = 0, target_found = 0;
    long long n_bases = 0;
    {
        const char *name, *seq;
        int32_t len;
        int r;
        while ((r = seqrec_next(rr, &name, &seq, &len)) > 0) {
            if (target_chr && strcmp(name, target_chr) != 0) continue;
            target_found = 1;
            contig_count++;
            n_bases += len;
            local.n_hit_tok = 0;
            scan_read_for_tokens(seq, len, d->key2i, local.pv, &local.hit_tok, &local.n_hit_tok, &local.m_hit_tok);
            if (local.n_hit_tok > 0) {
                int32_t cidx = add_contig(&contig_names, &n_contigs, &m_contigs, name);
                scan_contig(&shared, &local, hits, cidx, seq, len);
            }
            for (int32_t k = 0; k < local.n_hit_tok; ++k) local.pv[local.hit_tok[k]].n = 0;
            if (target_chr) break;
        }
        if (r < 0) {
            const char *e = seqsrc_error(src);
            fprintf(stderr, "[ERROR] %s: %s\n", o->input, *e ? e : "malformed FASTA/FASTQ record");
            rc = 2;
        }
    }
    seqrec_close(rr);
    if (seqsrc_close(src) != 0) rc = 2;
    if (rc) return rc;
    if (target_chr && !target_found) {
        fprintf(stderr, "[ERROR] sequence '%s' not found in %s\n", target_chr, o->input);
        return 2;
    }
    fprintf(stderr, "[INFO] processed %d record%s (%lld bases); %d hit the catalog\n",
            contig_count, contig_count == 1 ? "" : "s", n_bases, n_contigs);

    FILE *out = 0, *invert_out = 0;
    if (open_outputs(o, &out, &invert_out) != 0) return 2;

    FILE *err_out = 0;
    if (o->error_path && o->error_path[0]) {
        err_out = fopen(o->error_path, "wb");
        if (!err_out) { fprintf(stderr, "[ERROR] cannot open --many-kmer-hits file %s\n", o->error_path); return 2; }
        fprintf(err_out, "chr\tregion_start\tregion_end\tkmerA_id\tkmerB_id\t"
                         "kmerA_contig\tkmerA_location\tkmerA_hit_count\t"
                         "kmerB_contig\tkmerB_location\tkmerB_hit_count\n");
    }

    // -d: one row per unique pair saying whether it reached the main output
    // and, if not, why.
    gzFile dbg_out = 0;
    if (open_debug(o, d, &dbg_out) != 0) return 2;
    u64set_t *dbg_seen = dbg_out ? u64set_init() : 0;

    // pairs that did not resolve to one locus go to -m (--many-kmer-hits).
    if (err_out) {
        for (int32_t pi = 0; pi < n_pairs; ++pi) {
            pair_hit_info_t *info = &hits[pi];
            if (pair_hits_problem(info)) {
                /* the first contig each k-mer hit, if any */
                int32_t ca = info->A.contig_idx[0], cb = info->B.contig_idx[0];
                const char *a_contig = (ca >= 0 && ca < n_contigs) ? contig_names[ca] : "";
                const char *b_contig = (cb >= 0 && cb < n_contigs) ? contig_names[cb] : "";
                fprintf(err_out, "%s\t%lld\t%lld\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n",
                        pairs[pi].chr ? pairs[pi].chr : "",
                        (long long)pairs[pi].region_start, (long long)pairs[pi].region_end,
                        tok[info->A.idx].id, tok[info->B.idx].id,
                        a_contig, info->A.location[0], info->A.count,
                        b_contig, info->B.location[0], info->B.count);
            }
        }
    }

    for (int32_t pi = 0; pi < n_pairs; ++pi) {
        int32_t t1 = pairs[pi].tok1, t2 = pairs[pi].tok2;
        pair_hit_info_t *info = &hits[pi];
        const char *reason = pair_hits_problem(info);

        int n_emitted = 0;
        for (int s = 0; !reason && s < info->A.n_contigs; ++s) {
            const int sb = pair_hits_partner(info, s);
            const char *why = 0;
            int32_t posA = info->A.location[s], posB = sb >= 0 ? info->B.location[sb] : -1;
            if (sb < 0 || posA < 0 || posB < 0) why = "kmer location unavailable";

            if (!why) {
                int32_t left = posA < posB ? posA : posB;
                int32_t right = posA < posB ? posB : posA;
                if (right - (left + 31) < 0) why = "kmers overlap";
            }

            const anchor_slice_t *sl = &info->slice[s];
            if (!why && !sl->seq) why = "contig sequence unavailable";

            if (!why) {
                int32_t contig_idx = info->A.contig_idx[s];
                const char *emit_contig = (contig_idx >= 0 && contig_idx < n_contigs) ? contig_names[contig_idx] : "";
                char kmerA_strand = tok[t1].strand, kmerB_strand = tok[t2].strand;
                char kmerA_match = info->A.match_strand[s] ? info->A.match_strand[s] : tok[t1].strand;
                char kmerB_match = info->B.match_strand[sb] ? info->B.match_strand[sb] : tok[t2].strand;

                int xor_mismatch = (kmerA_strand != kmerA_match) != (kmerB_strand != kmerB_match);
                int both_match = (kmerA_strand == kmerA_match && kmerB_strand == kmerB_match);
                int both_mismatch = (kmerA_strand != kmerA_match && kmerB_strand != kmerB_match);
                int to_invert = xor_mismatch
                    || (both_match && posA > posB)
                    || (both_mismatch && posB > posA);

                if (!to_invert || invert_out)
                    record_hit(&pairs[pi], pi, to_invert ? 1 : 0, o->save_qnames ? emit_contig : 0, tok, t1, t2,
                               kmerA_match, kmerB_match, sl->posA, sl->posB, sl->seq, o->save_seq, 0, 0);

                if (to_invert) why = "strand/order mismatch";
                else n_emitted++;
            }
            if (why && !n_emitted && s == info->A.n_contigs - 1) reason = why;
        }

        if (dbg_out) {
            int absent;
            u64set_put(dbg_seen, pack_tt(t1, t2), &absent);
            if (absent)
                gzprintf(dbg_out, "%s\t%s\t%s\t%s\n", tok[t1].id, tok[t2].id,
                         reason ? "yes" : "no", reason ? reason : "NA");
        }
    }

    if (finish_outputs(o, pairs, n_pairs, tok, out, invert_out) != 0) rc = 2;
    if (err_out) fclose(err_out);
    if (dbg_out) gzclose(dbg_out);
    if (dbg_seen) u64set_destroy(dbg_seen);

    for (int32_t pi = 0; pi < n_pairs; ++pi)
        for (int s = 0; s < SEQ_MAX_CONTIGS_CAP; ++s) free(hits[pi].slice[s].seq);
    free(hits);
    for (int32_t i = 0; i < n_contigs; ++i) free(contig_names[i]);
    free(contig_names);
    scan_local_free(&local, &shared);
    scan_shared_free(&shared);
    return rc;
}

static int run_with_catalog(const opt_t *o, const char *chr_filter)
{
    pair_t *pairs = 0; int32_t n_pairs = 0, m_pairs = 0;
    token_t *tok = 0;  int32_t n_tok = 0, m_tok = 0;
    skipped_row_t *skip = 0; int32_t n_skip = 0, m_skip = 0;

    str2i_t *token_index = str2i_init();   /* "id\tstrand" -> tok idx */

    if (!path_exists(o->catalog)) {
        fprintf(stderr, "[ERROR] cannot find the k-mer catalog %s\n", o->catalog);
        return 2;
    }
    int rc_cat = read_pairs_file(o->catalog, chr_filter, &pairs, &n_pairs, &m_pairs,
                                 &tok, &n_tok, &m_tok, token_index,
                                 &skip, &n_skip, &m_skip);
    if (rc_cat != 0) {
        if (rc_cat == -1) fprintf(stderr, "[ERROR] cannot read k-mer catalog %s\n", o->catalog);
        return 2;   /* -2: the header problem was already reported */
    }

    if (n_pairs == 0) {
        if (chr_filter)
            fprintf(stderr, "[WARNING] no catalog pairs for %s in %s\n", chr_filter, o->catalog);
        else
            fprintf(stderr, "[WARNING] no pairs in catalog %s\n", o->catalog);
        return 0;
    }

    add_strand_twins(&tok, &n_tok, &m_tok, token_index);

    /* 31-mer key -> token index */
    key2i_t *key2i = key2i_init();

    if (build_dictionary(tok, n_tok, key2i) != 0) {
        fprintf(stderr, "[ERROR] failed building k-mer dictionary\n");
        return 2;
    }

    dict_t d;
    d.pairs = pairs; d.n_pairs = n_pairs;
    d.tok = tok;     d.n_tok = n_tok;
    d.key2i = key2i;
    d.token_index = token_index;
    d.skip = skip;   d.n_skip = n_skip;

    int rc = o->is_reads                ? run_reads(&d, o)
           : sr_fmt_is_alignment(o->fmt) ? run_alignment(&d, o)
           :                               run_sequence(&d, o);

    key2i_destroy(key2i);

    for (int32_t i = 0; i < n_tok; ++i) free(tok[i].id);
    free(tok);
    for (int32_t i = 0; i < n_pairs; ++i) {
        free(pairs[i].chr);
        free(pairs[i].loci);
        pair_free_extras(&pairs[i]);
    }
    free(pairs);

    for (int32_t i = 0; i < n_skip; ++i) {
        free(skip[i].idA);
        free(skip[i].idB);
    }
    free(skip);

    for (khint_t it = 0; it != kh_end(token_index); ++it)
        if (kh_exist(token_index, it)) free((char*)kh_key(token_index, it));
    str2i_destroy(token_index);

    return rc;
}

/* whole-genome runs, split by chromosome */

typedef struct {
    const char *chr;            /* chromosome as the catalog spells it */
    char  region[256];          /* reference name as this file's header spells it */
    char  out[PATH_MAX], err[PATH_MAX], inv[PATH_MAX], dbg[PATH_MAX];
    pid_t pid;
    int   ok;
} shard_t;

static char **samtools_header_refs(const char *path, const char *ref, int *n_out)
{
    int fds[2];
    *n_out = 0;
    if (pipe(fds) != 0) return 0;

    pid_t pid = fork();
    if (pid < 0) { close(fds[0]); close(fds[1]); return 0; }
    if (pid == 0) {
        dup2(fds[1], STDOUT_FILENO);
        close(fds[0]); close(fds[1]);
        const char *argv[8]; int n = 0;
        argv[n++] = "samtools"; argv[n++] = "view"; argv[n++] = "-H";
        if (ref && *ref) { argv[n++] = "-T"; argv[n++] = ref; }
        argv[n++] = path;
        argv[n] = 0;
        execvp("samtools", (char *const *)argv);
        _exit(127);
    }
    close(fds[1]);

    FILE *f = fdopen(fds[0], "r");
    char **names = 0; int n = 0, m = 0;
    if (f) {
        char line[8192];
        while (fgets(line, sizeof line, f)) {
            if (strncmp(line, "@SQ", 3) != 0) continue;
            char *sn = strstr(line, "\tSN:");
            if (!sn) continue;
            sn += 4;
            char *end = strpbrk(sn, "\t\r\n");
            if (end) *end = 0;
            if (n == m) { m = m ? m * 2 : 16; names = (char **)realloc(names, (size_t)m * sizeof *names); }
            names[n++] = strdup(sn);
        }
        fclose(f);
    } else close(fds[0]);

    int st = 0;
    while (waitpid(pid, &st, 0) < 0 && errno == EINTR) {}
    if (!WIFEXITED(st) || WEXITSTATUS(st) != 0) {
        for (int i = 0; i < n; ++i) free(names[i]);
        free(names);
        return 0;
    }
    *n_out = n;
    return names;
}

static int url_index_usable(const char *url, const char *ref, const char *refname)
{
    char region[300];
    snprintf(region, sizeof region, "%s:1-2", refname);

    fflush(NULL);
    pid_t pid = fork();
    if (pid < 0) return 0;
    if (pid == 0) {
        int devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) {
            dup2(devnull, STDOUT_FILENO);
            dup2(devnull, STDERR_FILENO);
            close(devnull);
        }
        const char *argv[8]; int n = 0;
        argv[n++] = "samtools"; argv[n++] = "view"; argv[n++] = "-c";
        if (ref && *ref) { argv[n++] = "-T"; argv[n++] = ref; }
        argv[n++] = url;
        argv[n++] = region;
        argv[n] = 0;
        execvp("samtools", (char *const *)argv);
        _exit(127);
    }

    int st = 0;
    while (waitpid(pid, &st, 0) < 0 && errno == EINTR) {}
    return WIFEXITED(st) && WEXITSTATUS(st) == 0;
}
static int url_is_shardable(const opt_t *o)
{
    int n = 0;
    char **refs = samtools_header_refs(o->input, o->ref, &n);
    if (!refs || n == 0) {
        for (int i = 0; i < n; ++i) free(refs[i]);
        free(refs);
        return -1;
    }
    int ok = url_index_usable(o->input, o->ref, refs[0]);
    for (int i = 0; i < n; ++i) free(refs[i]);
    free(refs);
    return ok;
}

/* fill in which chromosomes to run and what this file calls them. */
static int plan_shards(const opt_t *o, char *const *chrs, int n_chr, shard_t *sh)
{
    char **hdr_refs = 0; int n_hdr = 0;
    sr_t *rd = 0;

    if (o->fmt == SR_FMT_CRAM || o->fmt == SR_FMT_URL) {
        hdr_refs = samtools_header_refs(o->input, o->ref, &n_hdr);
        if (!hdr_refs)
            fprintf(stderr, "[INFO] cannot read the %s header with samtools; assuming catalog chromosome names\n", sr_fmt_name(o->fmt));
    } else {
        rd = sr_open(o->input, o->ref, 0, o->n_threads);
        if (!rd) { fprintf(stderr, "[ERROR] %s\n", sr_open_error()); return -1; }
    }

    int n = 0, skipped = 0;
    for (int i = 0; i < n_chr; ++i) {
        const char *name = chrs[i];

        char flip[256];
        if (strncasecmp(name, "chr", 3) == 0) snprintf(flip, sizeof flip, "%s", name + 3);
        else                                  snprintf(flip, sizeof flip, "chr%s", name);

        const char *pick = 0;
        if (hdr_refs) {
            for (int j = 0; j < n_hdr && !pick; ++j)
                if (chr_name_eq(hdr_refs[j], name)) pick = hdr_refs[j];
        } else if (rd) {
            if      (sr_has_ref(rd, name) == 1) pick = name;
            else if (sr_has_ref(rd, flip) == 1) pick = flip;
            else if (sr_has_ref(rd, name) < 0)  pick = name;
        } else {
            pick = name;
        }

        if (!pick) { ++skipped; continue; }
        sh[n].chr = name;
        snprintf(sh[n].region, sizeof sh[n].region, "%s", pick);
        sh[n].pid = -1;
        sh[n].ok  = 0;
        ++n;
    }

    if (rd) sr_close(rd);
    for (int i = 0; i < n_hdr; ++i) free(hdr_refs[i]);
    free(hdr_refs);

    if (skipped)
        fprintf(stderr, "[INFO] %d of %d catalog chromosomes are not in the %s header; skipped\n",
                skipped, n_chr, sr_fmt_name(o->fmt));
    return n;
}

static int try_mkdtemp(const char *pattern, char *out, size_t n)
{
    char tmpl[PATH_MAX];
    snprintf(tmpl, sizeof tmpl, "%s", pattern);
    if (!mkdtemp(tmpl)) return -1;
    snprintf(out, n, "%s", tmpl);
    return 0;
}

static int make_shard_dir(const opt_t *o, char *out, size_t n)
{
    char beside[PATH_MAX];
    beside[0] = 0;

    if (o->out_path && strcmp(o->out_path, "-") != 0) {
        const char *slash = strrchr(o->out_path, '/');
        if (slash) snprintf(beside, sizeof beside, "%.*s/.trumut.XXXXXX", (int)(slash - o->out_path), o->out_path);
        else       snprintf(beside, sizeof beside, ".trumut.XXXXXX");
        if (try_mkdtemp(beside, out, n) == 0) return 0;
    }

    const char *td = getenv("TMPDIR");
    if (!td || !*td) td = "/tmp";
    size_t tl = strlen(td);
    char fallback[PATH_MAX];
    snprintf(fallback, sizeof fallback, "%s%strumut.XXXXXX", td, (tl && td[tl-1] == '/') ? "" : "/");
    if (try_mkdtemp(fallback, out, n) == 0) return 0;

    fprintf(stderr, "[ERROR] cannot create a temporary directory (tried %s%s%s): %s\n",
            beside[0] ? beside : "", beside[0] ? " and " : "", fallback, strerror(errno));
    return -1;
}

static shard_t *g_sh = 0;
static int      g_n_sh = 0;
static char     g_dir[PATH_MAX];

static void shard_cleanup(void)
{
    if (!g_dir[0]) return;
    for (int i = 0; i < g_n_sh; ++i) {
        unlink(g_sh[i].out);
        unlink(g_sh[i].err);
        unlink(g_sh[i].inv);
        unlink(g_sh[i].dbg);
    }
    rmdir(g_dir);
    g_dir[0] = 0;
}

static void shard_signal(int sig)
{
    for (int i = 0; i < g_n_sh; ++i)
        if (g_sh[i].pid > 0) kill(g_sh[i].pid, SIGTERM);

    shard_cleanup();
    signal(sig, SIG_DFL);
    raise(sig);           /* die of the original signal, so callers see it */
}

static int append_file(FILE *dst, const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    char buf[1 << 16];
    size_t n;
    int rc = 0;
    while ((n = fread(buf, 1, sizeof buf, f)) > 0)
        if (fwrite(buf, 1, n, dst) != n) { rc = -1; break; }
    fclose(f);
    return rc;
}

/* header once, then every chromosome's rows in karyotype order. */
static int merge_one(const char *dest, const char *header, shard_t *sh, int n_sh,
                     size_t off, const char *what)
{
    FILE *out = stdout;
    if (dest && strcmp(dest, "-") != 0) {
        out = fopen(dest, "wb");
        if (!out) { fprintf(stderr, "[ERROR] cannot open %s %s\n", what, dest); return -1; }
    }
    if (header) fputs(header, out);

    int rc = 0;
    for (int i = 0; i < n_sh; ++i)
        if (sh[i].ok && append_file(out, (const char *)&sh[i] + off) != 0)
            rc = -1;

    if (out != stdout) fclose(out);
    else fflush(out);
    if (rc != 0) fprintf(stderr, "[ERROR] writing %s failed\n", what);
    return rc;
}

static int merge_one_gz(const char *dest, const char *header, shard_t *sh, int n_sh,
                        size_t off, const char *what)
{
    FILE *out;
    if (dest && strcmp(dest, "-") != 0) {
        out = fopen(dest, "wb");
        if (!out) { fprintf(stderr, "[ERROR] cannot open %s %s\n", what, dest); return -1; }
    } else {
        out = stdout;
    }

    int rc = 0;
    fflush(out);
    gzFile g = gzdopen(dup(fileno(out)), "wb");
    if (!g || gzputs(g, header) < 0) rc = -1;
    if (g) gzclose(g);

    for (int i = 0; i < n_sh && rc == 0; ++i)
        if (sh[i].ok && append_file(out, (const char *)&sh[i] + off) != 0)
            rc = -1;

    if (out != stdout) fclose(out);
    else fflush(out);
    if (rc != 0) fprintf(stderr, "[ERROR] writing %s failed\n", what);
    return rc;
}

static void free_chr_list(char **chrs, int n)
{
    for (int i = 0; i < n; ++i) free(chrs[i]);
    free(chrs);
}

static int chr_ptr_cmp(const void *a, const void *b)
{
    return chr_order_cmp(*(char *const *)a, *(char *const *)b);
}

static int run_sharded(const opt_t *o)
{
    char **chrs = 0;
    int n_chr = catalog_chr_list(o->catalog, &chrs);
    if (n_chr < 0) {
        fprintf(stderr, "[ERROR] cannot read k-mer catalog %s\n", o->catalog);
        return 2;
    }
    if (n_chr == 0) {
        fprintf(stderr, "[ERROR] k-mer catalog %s has no rows\n", o->catalog);
        free_chr_list(chrs, n_chr);
        return 2;
    }

    qsort(chrs, (size_t)n_chr, sizeof *chrs, chr_ptr_cmp);

    shard_t *sh = (shard_t *)calloc(n_chr, sizeof *sh);
    if (!sh) { fprintf(stderr, "[ERROR] out of memory\n"); free_chr_list(chrs, n_chr); return 2; }

    int n_sh = plan_shards(o, chrs, n_chr, sh);
    if (n_sh < 0)  { free(sh); free_chr_list(chrs, n_chr); return 2; }
    if (n_sh == 0) {
        fprintf(stderr, "[ERROR] none of the catalog's %d chromosomes are in the %s header; check the reference names or pass -t\n", n_chr, sr_fmt_name(o->fmt));
        free(sh);
        free_chr_list(chrs, n_chr);
        return 2;
    }

    char why[256];
    long long per_worker = worker_mem_estimate(o->catalog, n_sh);
    int n_jobs = o->n_jobs > 0 ? o->n_jobs : auto_jobs(per_worker, n_sh, why, sizeof why);
    if (n_jobs > n_sh) n_jobs = n_sh;

    const char *plural = n_sh == 1 ? "chromosome" : "chromosomes";
    if (n_jobs == 1)
        fprintf(stderr, "[INFO] running %d %s, one at a time%s\n",
                n_sh, plural, o->n_jobs > 0 ? " (-j 1)" : "");
    else if (o->n_jobs > 0)
        fprintf(stderr, "[INFO] running %d %s, %d at a time (-j %d)\n", n_sh, plural, n_jobs, n_jobs);
    else
        fprintf(stderr, "[INFO] running %d %s, %d at a time (%s); use -j to change\n",
                n_sh, plural, n_jobs, why);

    int per_worker_threads = o->n_threads;

    char dir[PATH_MAX];
    if (make_shard_dir(o, dir, sizeof dir) != 0) { free(sh); free_chr_list(chrs, n_chr); return 2; }

    for (int i = 0; i < n_sh; ++i) {
        const char *ogz = wants_gz(o->out_path)    ? ".gz" : "";
        const char *igz = wants_gz(o->invert_path) ? ".gz" : "";
        if (snprintf(sh[i].out, sizeof sh[i].out, "%s/%s.tsv%s",     dir, sh[i].region, ogz) >= (int)sizeof sh[i].out ||
            snprintf(sh[i].err, sizeof sh[i].err, "%s/%s.err.tsv",   dir, sh[i].region) >= (int)sizeof sh[i].err ||
            snprintf(sh[i].inv, sizeof sh[i].inv, "%s/%s.inv.tsv%s", dir, sh[i].region, igz) >= (int)sizeof sh[i].inv ||
            snprintf(sh[i].dbg, sizeof sh[i].dbg, "%s/%s.dbg.tsv", dir, sh[i].region) >= (int)sizeof sh[i].dbg) {
            fprintf(stderr, "[ERROR] the temporary path for %s does not fit in %d bytes: %s\n",
                    sh[i].region, PATH_MAX, dir);
            rmdir(dir);
            free(sh);
            free_chr_list(chrs, n_chr);
            return 2;
        }
    }

    int want_err = o->error_path  && o->error_path[0];
    int want_inv = o->invert_path && o->invert_path[0];
    int want_dbg = o->debug_path  && o->debug_path[0];

    g_sh = sh; g_n_sh = n_sh;
    snprintf(g_dir, sizeof g_dir, "%s", dir);
    struct sigaction sa, old_int, old_term, old_hup, old_pipe;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = shard_signal;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGINT,  &sa, &old_int);
    sigaction(SIGTERM, &sa, &old_term);
    sigaction(SIGHUP,  &sa, &old_hup);
    sigaction(SIGPIPE, &sa, &old_pipe);

    int next = 0, running = 0, n_failed = 0, n_done = 0;

    while (next < n_sh || running > 0) {
        while (running < n_jobs && next < n_sh) {
            fflush(NULL);
            pid_t pid = fork();
            if (pid < 0) {
                fprintf(stderr, "[ERROR] cannot fork a worker for %s: %s\n",
                        sh[next].region, strerror(errno));
                sh[next].pid = -1;
                ++next; ++n_failed;
                continue;
            }
            if (pid == 0) {
                g_dir[0] = 0;
                signal(SIGINT, SIG_DFL); signal(SIGTERM, SIG_DFL);
                signal(SIGHUP, SIG_DFL); signal(SIGPIPE, SIG_DFL);

                opt_t co = *o;
                co.target_chr  = sh[next].region;
                co.out_path    = sh[next].out;
                co.error_path  = want_err ? sh[next].err : 0;
                co.invert_path = want_inv ? sh[next].inv : 0;
                co.debug_path  = want_dbg ? sh[next].dbg : 0;
                co.n_threads   = per_worker_threads;
                co.shard       = 1;
                _exit(run_with_catalog(&co, sh[next].chr));
            }
            sh[next].pid = pid;
            ++running; ++next;
        }

        int st = 0;
        pid_t done = waitpid(-1, &st, 0);
        if (done < 0) {
            if (errno == EINTR) continue;
            break;
        }
        for (int i = 0; i < n_sh; ++i) {
            if (sh[i].pid != done) continue;
            sh[i].pid = -1;
            --running;
            ++n_done;
            int rc = WIFEXITED(st) ? WEXITSTATUS(st) : -1;
            if (rc == 0) {
                sh[i].ok = 1;
                fprintf(stderr, "  [%2d/%d] %s done\n", n_done, n_sh, sh[i].region);
            } else {
                ++n_failed;
                if (WIFSIGNALED(st))
                    fprintf(stderr, "  [%2d/%d] %s FAILED (killed by signal %d)\n",
                            n_done, n_sh, sh[i].region, WTERMSIG(st));
                else
                    fprintf(stderr, "  [%2d/%d] %s FAILED (exit %d)\n",
                            n_done, n_sh, sh[i].region, rc);
            }
            break;
        }
    }

    /* merge what succeeded, then report failures. */
    int rc = 0;
    int (*merge_out)(const char *, const char *, shard_t *, int, size_t, const char *) =
        wants_gz(o->out_path) ? merge_one_gz : merge_one;
    int (*merge_inv)(const char *, const char *, shard_t *, int, size_t, const char *) =
        wants_gz(o->invert_path) ? merge_one_gz : merge_one;
    if (merge_out(o->out_path, header_for(o), sh, n_sh, offsetof(shard_t, out), "output") != 0) rc = 2;
    if (want_inv && merge_inv(o->invert_path, header_for(o), sh, n_sh, offsetof(shard_t, inv), "invert file") != 0) rc = 2;
    if (want_err && merge_one(o->error_path, DETECT_ERROR_HEADER, sh, n_sh, offsetof(shard_t, err), "many-kmer-hits file") != 0) rc = 2;
    if (want_dbg && merge_one_gz(o->debug_path, DETECT_DEBUG_HEADER, sh, n_sh, offsetof(shard_t, dbg), "debug file") != 0) rc = 2;

    shard_cleanup();
    sigaction(SIGINT,  &old_int,  0);
    sigaction(SIGTERM, &old_term, 0);
    sigaction(SIGHUP,  &old_hup,  0);
    sigaction(SIGPIPE, &old_pipe, 0);
    g_sh = 0; g_n_sh = 0;

    if (n_failed) {
        fprintf(stderr, "[ERROR] %d of %d chromosomes failed; output is incomplete (%d succeeded), see messages above\n",
                n_failed, n_sh, n_sh - n_failed);
        rc = 2;
    }

    free(sh);
    free_chr_list(chrs, n_chr);
    return rc;
}

static void allele_usage(FILE *fp)
{
    fprintf(fp,
    "Usage: trumut allele <input> <catalog> <output> --input-type reads|assembly [options]\n"
    "\n"
    "Positional:\n"
    "  <input>     alignments (SAM, SAM.gz, BAM, CRAM), an assembly (FASTA), or unaligned\n"
    "              reads (FASTQ, FASTA, unaligned SAM/BAM); a local file, stdin (-), or a URL\n"
    "  <catalog>   uniqmer-annotated tandem repeat catalog (pre-built, or from trumut prep-catalog)\n"
    "  <output>    output path\n"
    "\n"
    "Required:\n"
    "  --input-type STR      reads for sequencing reads, assembly for assembled contigs\n"
    "  --ploidy N            required with assembly: haplotype contigs per locus (1 haploid, 2 diploid)\n"
    "\n"
    "Input:\n"
    "  --input-fmt STR       specify input file format [sniffed]\n"
    "  -r, --reference FILE  reference FASTA, needed for CRAM\n"
    "\n"
    "Output:\n"
    "  --save-allele-seq     add an allele_seqs column: the read sequence between the two uniqmers\n"
    "  --save-query-names    add a query_names column: read names, or contig names for assembly input\n"
    "  --save-kmer-coords    add kmerA_ref_coords and kmerB_ref_coords columns: the reference position of each uniqmer\n"
    "  -w, --wrong-kmer-hits FILE  uniqmer pairs found in a read / contig on unexpected strands or in the wrong order\n"
    "  -m, --many-kmer-hits FILE   uniqmer pairs found more than once in a single read / contig\n"
    "  -d, --debug-file FILE       one row per catalog uniqmer pair: whether it was dropped and, if so, why\n"
    "\n"
    "Indexed alignments (BAM, CRAM, or URL with a remote index):\n"
    "  -t, --target STR      process one chromosome only, spelled exactly as in the BAM/CRAM header [all]\n"
    "  -j, --jobs INT        chromosomes processed in parallel, one worker each; lowered if they would not fit in memory [3]\n"
    "  -q, --min-mapq INT    skip alignments below this MAPQ [0]\n"
    "  -@, --threads INT     threads passed to samtools for CRAM decoding, per worker [3]\n"
    "\n"
    "Unindexed streams (SAM, unindexed BAM, stdin, FASTQ/FASTA reads):\n"
    "  -t, --target STR      SAM/BAM streams only: process one chromosome; chr1 and 1 both match [all]\n"
    "  -j, --jobs INT        scanning threads [3]\n"
    "  -q, --min-mapq INT    SAM/BAM streams only: skip alignments below this MAPQ [0]\n"
    "  -@, --threads INT     BGZF decompression threads for bgzip FASTQ or unaligned BAM [3]\n"
    "\n"
    "Assemblies (FASTA contigs):\n"
    "  -t, --target STR      process one record (contig) only, by name [all]\n"
    "  -@, --threads INT     BGZF decompression threads for a bgzip assembly [3]\n"
    "\n"
    "  -h, --help            show this message\n"
    "\n"
    "Documentation: https://erkim11.github.io/trumut/\n");
}

static void reads_thread_plan(opt_t *o)
{
    o->inflate_threads = o->n_threads;
    if (o->inflate_threads < 1) o->inflate_threads = 1;
}

int allele_main(int argc, char *argv[])
{
    opt_t o;
    memset(&o, 0, sizeof o);
    o.n_threads = 3;   /* -@ default */
    o.fmt = SR_FMT_UNKNOWN;

    const char *fmt_arg = 0;
    int saw_ref = 0, saw_mapq = 0, saw_jobs = 0, saw_ploidy = 0;

    static ko_longopt_t longopts[] = {
        { (char*)"input-fmt", ko_required_argument, 301 },
        { (char*)"save-query-names", ko_no_argument, 302 },
        { (char*)"save-allele-seq",  ko_no_argument, 303 },
        { (char*)"save-kmer-coords", ko_no_argument, 305 },
        { (char*)"ploidy",           ko_required_argument, 304 },
        { (char*)"input-type",       ko_required_argument, 306 },
        { (char*)"target",    ko_required_argument, 't' },
        { (char*)"many-kmer-hits", ko_required_argument, 'm' },
        { (char*)"wrong-kmer-hits", ko_required_argument, 'w' },
        { (char*)"jobs",      ko_required_argument, 'j' },
        { (char*)"reference", ko_required_argument, 'r' },
        { (char*)"min-mapq",  ko_required_argument, 'q' },
        { (char*)"threads",   ko_required_argument, '@' },
        { (char*)"debug-file", ko_required_argument, 'd' },
        { (char*)"help",      ko_no_argument,       'h' },
        { 0, 0, 0 }
    };

    ketopt_t k = KETOPT_INIT;
    int c;
    while ((c = ketopt(&k, argc, argv, 1, "r:m:t:w:q:@:j:d:h", longopts)) >= 0) {
        if      (c == 'r') { o.ref = k.arg; saw_ref = 1; }
        else if (c == 't') o.target_chr = k.arg;
        else if (c == 'm') o.error_path = k.arg;   /* --many-kmer-hits */
        else if (c == 'w') o.invert_path = k.arg;
        else if (c == 'q') { if (parse_int_opt("-q", k.arg, &o.mapq_min) != 0) return 2; saw_mapq = 1; }
        else if (c == '@') { if (parse_int_opt("-@", k.arg, &o.n_threads) != 0) return 2; }
        else if (c == 'j') {
            if (parse_int_opt("-j", k.arg, &o.n_jobs) != 0) return 2;
            saw_jobs = 1;
            if (o.n_jobs < 1) {
                fprintf(stderr, "[ERROR] -j must be at least 1 (got '%s')\n", k.arg);
                return 2;
            }
        }
        else if (c == 'd') o.debug_path = k.arg;
        else if (c == 301) fmt_arg = k.arg;
        else if (c == 302) o.save_qnames = 1;
        else if (c == 303) o.save_seq = 1;
        else if (c == 305) o.save_coords = 1;
        else if (c == 306) {
            if      (strcmp(k.arg, "reads") == 0)    o.reads_mode = 1;
            else if (strcmp(k.arg, "assembly") == 0) o.reads_mode = -1;
            else {
                fprintf(stderr, "[ERROR] --input-type must be 'reads' or 'assembly' (got '%s')\n", k.arg);
                return 2;
            }
        }
        else if (c == 304) {
            if (parse_int_opt("--ploidy", k.arg, &o.ploidy) != 0) return 2;
            saw_ploidy = 1;
            if (o.ploidy < 1 || o.ploidy > SEQ_MAX_CONTIGS_CAP) {
                fprintf(stderr, "[ERROR] --ploidy must be between 1 and %d (got '%s')\n", SEQ_MAX_CONTIGS_CAP, k.arg);
                return 2;
            }
        }
        else if (c == 'h') { allele_usage(stdout); return 0; }
        else if (c == '?') { fprintf(stderr, "[ERROR] unknown option\n"); allele_usage(stderr); return 2; }
        else if (c == ':') { fprintf(stderr, "[ERROR] missing argument\n"); allele_usage(stderr); return 2; }
    }

    seq_max_contigs = o.ploidy;

    /* three positionals, in order: <input> <catalog> <output>. */
    const int n_pos = argc - k.ind;
    if (n_pos != 3) {
        if (n_pos == 0) { allele_usage(stderr); return 2; }
        if (n_pos < 3)
            fprintf(stderr, "[ERROR] allele takes three arguments: <input> <catalog> <output> (got %d)\n", n_pos);
        else
            fprintf(stderr, "[ERROR] unexpected extra argument '%s'; allele takes <input> <catalog> <output>\n", argv[k.ind + 3]);
        allele_usage(stderr);
        return 2;
    }
    o.input    = argv[k.ind];
    o.catalog  = argv[k.ind + 1];
    o.out_path = argv[k.ind + 2];

    /* --input-type: reads vs assembly */
    if (o.reads_mode == 0) {
        fprintf(stderr, "[ERROR] --input-type is required: reads or assembly\n");
        return 2;
    }
    if (o.reads_mode == -1 && !saw_ploidy) {
        fprintf(stderr, "[ERROR] --ploidy is required with --input-type assembly (1 haploid, 2 diploid, ...)\n");
        return 2;
    }
    if (saw_ploidy && o.reads_mode != -1) {
        fprintf(stderr, "[ERROR] --ploidy applies to --input-type assembly only\n");
        return 2;
    }

    const int is_url = sr_path_is_url(o.input);
    const int is_stdin = strcmp(o.input, "-") == 0;
    reads_thread_plan(&o);

    if (fmt_arg) {
        o.fmt = sr_fmt_parse(fmt_arg);
        if (o.fmt == SR_FMT_UNKNOWN) {
            fprintf(stderr, "[ERROR] unknown --input-fmt '%s'; expected sam, bam, cram, fasta, or fastq\n", fmt_arg);
            return 2;
        }
    }
    const int fmt_says_seq = o.fmt == SR_FMT_FASTA || o.fmt == SR_FMT_FASTQ;
    int sniff = 0;
    if (is_stdin) sniff = 1;
    else if (is_url) {
        const int kind = seqsrc_url_kind(o.input);
        const int want_seq = o.reads_mode == -1 || fmt_says_seq;
        if (fmt_arg && !fmt_says_seq) o.fmt = SR_FMT_URL;
        else if (want_seq || kind == 1) sniff = 1;
        else o.fmt = SR_FMT_URL;
    } else if (!fmt_arg) {
        o.fmt = sr_probe_format(o.input);
        if (o.fmt == SR_FMT_UNKNOWN) {
            fprintf(stderr, "[ERROR] cannot determine the format of %s (unreadable or empty); pass --input-fmt\n", o.input);
            return 2;
        }
    }
    if (sniff) {
        char err[512];
        seqsrc_t *src = seqsrc_open(o.input, o.inflate_threads, err, sizeof err);
        if (!src) { fprintf(stderr, "[ERROR] %s\n", err); return 2; }
        uint8_t m[16];
        ssize_t n = seqsrc_peek(src, m, sizeof m);
        if (n <= 0) {
            fprintf(stderr, "[ERROR] %s: %s\n", is_stdin ? "stdin" : o.input,
                    n < 0 ? seqsrc_error(src) : "empty input");
            seqsrc_close(src);
            return 2;
        }
        if (!fmt_arg) o.fmt = sr_classify_bytes(m, (size_t)n);
        if (is_url && sr_fmt_is_alignment(o.fmt)) {
            seqsrc_close(src);
            o.fmt = SR_FMT_URL;
        } else {
            o.pre_src = src;
        }
    }
    const char *shown_input = is_stdin ? "stdin" : o.input;

    if (o.reads_mode == -1 && o.fmt != SR_FMT_FASTA) {
        if (o.fmt == SR_FMT_FASTQ)
            fprintf(stderr, "[ERROR] --input-type assembly needs a FASTA; %s is FASTQ (use --input-type reads, or convert with seqtk seq -A)\n", shown_input);
        else if (o.fmt == SR_FMT_URL)
            fprintf(stderr, "[ERROR] --input-type assembly needs a FASTA; %s holds alignments (use --input-type reads)\n", o.input);
        else
            fprintf(stderr, "[ERROR] --input-type assembly needs a FASTA; %s is %s (use --input-type reads)\n", shown_input, sr_fmt_name(o.fmt));
        if (o.pre_src) seqsrc_close(o.pre_src);
        return 2;
    }
    o.is_reads = o.reads_mode == 1 && (o.fmt == SR_FMT_FASTA || o.fmt == SR_FMT_FASTQ);
    if (o.is_reads || o.reads_mode == -1) {
        const char *what = o.is_reads ? "unaligned reads" : "an assembly";
        const char *bad = 0;
        if (o.is_reads && o.target_chr) {
            fprintf(stderr, "[ERROR] -t needs aligned input; unaligned reads carry no chromosome\n");
            if (o.pre_src) seqsrc_close(o.pre_src);
            return 2;
        }
        if (saw_ref) bad = "-r";
        else if (saw_mapq) bad = "-q";
        else if (saw_jobs && !o.is_reads) bad = "-j";
        if (bad) {
            fprintf(stderr, "[ERROR] %s applies to alignment input only, but %s is %s (%s)\n",
                    bad, shown_input, what, sr_fmt_name(o.fmt));
            if (o.pre_src) seqsrc_close(o.pre_src);
            return 2;
        }
    }
    if (o.is_reads) {
        fprintf(stderr, "[INFO] loading the whole catalog (unaligned reads carry no chromosome)\n");
        return run_with_catalog(&o, 0);
    }

    int is_aln = sr_fmt_is_alignment(o.fmt);

    /* a local BAM with no index (typically an unaligned BAM) cannot be split by chr, so it streams like SAM. */
    int bam_unindexed = 0;
    if (o.fmt == SR_FMT_BAM && !is_stdin) {
        char idx[PATH_MAX];
        bam_unindexed = find_index(o.input, ".bai", ".bam", idx, sizeof idx) != 0;
    }
    if (!bam_unindexed && require_index(&o) != 0) return 2;

    const char *chr_filter = (is_aln && o.target_chr) ? o.target_chr : 0;

    int is_stream = 0;
    if (is_aln && (o.fmt == SR_FMT_SAM || bam_unindexed || is_stdin)) is_stream = 1;
    if (bam_unindexed && !o.target_chr)
        fprintf(stderr, "[INFO] %s has no .bai index; streaming it once with the whole catalog (samtools index enables per-chromosome runs)\n", o.input);
    if (is_aln && o.fmt == SR_FMT_URL && !o.target_chr) {
        int shardable = url_is_shardable(&o);
        if (shardable == 1)
            return run_sharded(&o);
        if (shardable == 0)
            fprintf(stderr, "[INFO] %s has no remote index; streaming it once with the whole catalog (upload the .bai/.crai to split by chromosome)\n", o.input);
        is_stream = 1;
    }
    if (is_stream) {
        o.is_reads = 1;
        int n_scan = o.n_jobs > 0 ? o.n_jobs : DEFAULT_JOBS;
        if (is_stdin && !o.target_chr)
            fprintf(stderr, "[INFO] stdin (%s) cannot be split by chromosome; streaming once with the whole catalog on %d thread%s (-t limits the catalog)\n",
                sr_fmt_name(o.fmt), n_scan, n_scan == 1 ? "" : "s");
        else if (o.fmt == SR_FMT_SAM && !o.target_chr)
            fprintf(stderr, "[INFO] %s is SAM and cannot be split by chromosome; streaming once with the whole catalog on %d thread%s (an indexed BAM enables per-chromosome runs)\n",
                o.input, n_scan, n_scan == 1 ? "" : "s");
        else if (!o.target_chr)
            fprintf(stderr, "[INFO] stream input cannot be split by chromosome; streaming once with the whole catalog on %d thread%s (-t limits the catalog)\n",
                n_scan, n_scan == 1 ? "" : "s");
        return run_with_catalog(&o, chr_filter);
    }

    if (is_aln && !o.target_chr) {
        return run_sharded(&o);
    } else if (saw_jobs && is_aln) {
        fprintf(stderr, "[INFO] -j has no effect with -t\n");
    }

    if (!chr_filter && !is_aln && o.target_chr)
        fprintf(stderr, "[INFO] loading the whole catalog (-t names a %s record, not a chromosome)\n",
                sr_fmt_name(o.fmt));
    else if (!chr_filter)
        fprintf(stderr, "[INFO] loading the whole catalog\n");

    return run_with_catalog(&o, chr_filter);
}
