/*
catalog.h
*/

#ifndef CATALOG_H
#define CATALOG_H

#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "third_party/khashl.h"
#include "third_party/nt4-table.h"
#include "third_party/seqtk-revcomp.h"

#ifdef __cplusplus
#include <string>
#include <string_view>
#endif

KHASHL_MAP_INIT(static inline, str2i_t, str2i, const char*, int32_t, kh_hash_str, kh_eq_str)
KHASHL_SET_INIT(static inline, u64set_t, u64set, uint64_t, kh_hash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(static inline, key2i_t, key2i, uint64_t, int32_t, kh_hash_uint64, kh_eq_generic)

typedef struct {
    char   *id;       // original kmer_id
    char    strand;   // '+' or '-'
    char    seq[32];  // plus-strand 31-mer from the catalog, uppercased
    uint64_t key;     // 2-bit encoded 31-mer (revcomp applied when strand is '-')
} token_t;

// one distinct k-mer pair
typedef struct {
    int32_t tok1, tok2;   // indices into the token array
    int64_t pos1, pos2;   // hg38 positions (0-based start)
    char *chr;
    int32_t chr_rank;     // karyotype rank of chr (catalog.h), cached for sorting
    int64_t region_start; // min tr_start across the pair's loci
    int64_t region_end;   // max tr_end across the pair's loci
    char *loci;           // "start-end:motif" tokens, comma-separated
    void *motif_kmers;    // the pair's motif 7-mer set (std::unordered_set<uint64_t>*)
    void *acc[2];         // per-read values collected for the main output [0] and
                          // the -w invert output [1]
} pair_t;

void pair_free_extras(pair_t *p);

// hit positions of one token in one read.
typedef struct {
    int32_t n, m;
    int32_t *a; // positions in read
} posvec_t;

// the pairs a token belongs to.
typedef struct {
    int32_t n, m;
    int32_t *pairs;
} token_pairs_t;

// encode a 31-mer into a 2-bit key.
int encode_31mer(const char *seq31, uint64_t *out_key);

// split a line on tabs in place.
int split_tab(char *s, char **f, int max_f);

typedef struct {
    char *idA, *idB;
    const char *reason;
} skipped_row_t;

int chr_name_eq(const char *a, const char *b);

int read_pairs_file(const char *pairs_path, const char *chr_filter,
                    pair_t **pairs, int32_t *n_pairs, int32_t *m_pairs,
                    token_t **tok, int32_t *n_tok, int32_t *m_tok,
                    str2i_t *token_index,
                    skipped_row_t **skip, int32_t *n_skip, int32_t *m_skip);

// give every token its opposite-strand twin.
void add_strand_twins(token_t **tok, int32_t *n_tok, int32_t *m_tok,
                      str2i_t *token_index);

// encode every token's 31-mer and index it as key -> token.
int build_dictionary(token_t *tok, int32_t n_tok, key2i_t *key2i);

const char *detect_header(int with_qname, int with_coords, int with_seq);

// build every pair's motif 7-mer set up front.
void pair_precompute_motifs(pair_t *pairs, int32_t n_pairs);

// record one read's hit on a pair.
void record_hit(pair_t *pair, int32_t pair_index, int which,
                const char *name,
                const token_t *tok, int32_t t1_ref, int32_t t2_ref,
                char kmerA_match, char kmerB_match,
                int32_t posA, int32_t posB,
                const char *seq, int with_seq,
                const char *hp, const char *ps);

int flush_pair_rows(FILE *out, pair_t *pairs, int32_t n_pairs, const token_t *tok,
                    int which, int with_qname, int with_coords, int with_seq);

void scan_read_for_tokens(const char *seq, int32_t len, key2i_t *key2i,
                          posvec_t *pv,
                          int32_t **hit_tok, int32_t *n_hit_tok, int32_t *m_hit_tok);

static inline int32_t chr_rank_n(const char *name, size_t n)
{
    if (n >= 3 &&
        (name[0] == 'c' || name[0] == 'C') &&
        (name[1] == 'h' || name[1] == 'H') &&
        (name[2] == 'r' || name[2] == 'R')) { name += 3; n -= 3; }
    if (n == 0) return 1000;
    if (n <= 2) {
        int numeric = 1;
        int32_t v = 0;
        for (size_t i = 0; i < n; ++i) {
            char c = name[i];
            if (c < '0' || c > '9') { numeric = 0; break; }
            v = v * 10 + (c - '0');
        }
        if (numeric && v >= 1 && v <= 22) return v;
    }
    if (n == 1 && (name[0] == 'X' || name[0] == 'x')) return 23;
    if (n == 1 && (name[0] == 'Y' || name[0] == 'y')) return 24;
    if ((n == 1 && (name[0] == 'M' || name[0] == 'm')) ||
        (n == 2 && (name[0] == 'M' || name[0] == 'm') && (name[1] == 'T' || name[1] == 't')))
        return 25;
    return 1000;
}

static inline int32_t chr_rank_cstr(const char *name)
{
    return chr_rank_n(name, strlen(name));
}

static inline int chr_order_cmp(const char *a, const char *b)
{
    int32_t ra = chr_rank_cstr(a), rb = chr_rank_cstr(b);
    if (ra != rb) return ra < rb ? -1 : 1;
    return strcmp(a, b);
}

#ifdef __cplusplus
inline int32_t chr_rank(std::string_view name) { return chr_rank_n(name.data(), name.size()); }

struct ChrKey {
  int32_t rank = 0;
  std::string name;

  bool operator<(const ChrKey& o) const {
    if (rank != o.rank) return rank < o.rank;
    return name < o.name;
  }
  bool operator==(const ChrKey& o) const { return name == o.name; }
};
#endif

constexpr int MOTIF_KMER_K = 7;

// encode a k-mer string to its canonical 2-bit key.
inline uint64_t kmer_canonical_key(const char *seq, int k) {
  const uint64_t mask = (1ULL << (k * 2)) - 1;
  const int shift = (k - 1) * 2;
  uint64_t x0 = 0, x1 = 0;
  for (int i = 0; i < k; ++i) {
    int c = (int)kom_nt4_table[(uint8_t)seq[i]];
    if (c >= 4) return UINT64_MAX;
    x0 = (x0 << 2 | (uint64_t)c) & mask;
    x1 = (x1 >> 2) | ((uint64_t)(3 - c) << shift);
  }
  return x0 < x1 ? x0 : x1;
}

// canonical 7-mer set of one repeated motif.
inline std::unordered_set<uint64_t> find_motif_kmers(std::string_view motif) {
  const int k = MOTIF_KMER_K;
  std::unordered_set<uint64_t> out;
  if (motif.empty()) return out;
  const int m = (int)motif.size();
  out.reserve(32);
  std::string rep;
  if (m >= k) {
    rep.reserve(2 * m);
    rep.append(motif);
    rep.append(motif);
  } else {
    int nrep = (2 * k / m) + 1;
    rep.reserve(m * nrep);
    for (int i = 0; i < nrep; ++i) rep.append(motif);
  }
  const int n = (int)rep.size();
  for (int i = 0; i + k <= n; ++i) {
    uint64_t key = kmer_canonical_key(rep.data() + i, k);
    if (key != UINT64_MAX) out.insert(key);
  }
  return out;
}

inline std::unordered_set<uint64_t> find_loci_motif_kmers(std::string_view loci) {
  std::unordered_set<uint64_t> out;
  size_t pos = 0;
  while (pos < loci.size()) {
    size_t comma = loci.find(',', pos);
    if (comma == std::string_view::npos) comma = loci.size();
    size_t colon = loci.find(':', pos);
    if (colon != std::string_view::npos && colon < comma) {
      auto mk = find_motif_kmers(loci.substr(colon + 1, comma - colon - 1));
      out.insert(mk.begin(), mk.end());
    }
    pos = comma + 1;
  }
  return out;
}

inline uint64_t count_motif_kmers(const char *seq, int len,
                                  const std::unordered_set<uint64_t> &motifs) {
  const int k = MOTIF_KMER_K;
  const uint64_t mask = (1ULL << (k * 2)) - 1;
  const int shift = (k - 1) * 2;
  uint64_t x0 = 0, x1 = 0, hits = 0;
  int l = 0;
  for (int i = 0; i < len; ++i) {
    int c = (int)kom_nt4_table[(uint8_t)seq[i]];
    if (c < 4) {
      x0 = (x0 << 2 | (uint64_t)c) & mask;
      x1 = (x1 >> 2) | ((uint64_t)(3 - c) << shift);
      if (++l >= k && motifs.count(x0 < x1 ? x0 : x1)) hits++;
    } else { l = 0; x0 = x1 = 0; }
  }
  return hits;
}

inline double parse_motif_count(std::string_view sv) {
  if (sv.empty() || sv == "." || sv == "NA" || sv == "nan") return std::numeric_limits<double>::quiet_NaN();
  double v = 0; int frac = -1; size_t i = 0;
  for (; i < sv.size(); ++i) {
    char c = sv[i];
    if (c >= '0' && c <= '9') {
      if (frac < 0) v = v * 10 + (c - '0');
      else { v += (c - '0') / std::pow(10.0, ++frac); }
    } else if (c == '.' && frac < 0) {
      frac = 0;
    } else {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }
  return v;
}

inline double median_of(std::vector<double> &vals) {
  vals.erase(std::remove_if(vals.begin(), vals.end(),
                            [](double v) { return std::isnan(v); }),
             vals.end());
  if (vals.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(vals.begin(), vals.end());
  const size_t n = vals.size();
  return (n % 2) ? vals[n / 2] : (vals[n / 2 - 1] + vals[n / 2]) / 2.0;
}

inline std::string format_median(double v, const char *na = "NA") {
  if (std::isnan(v)) return na;
  char b[64];
  std::snprintf(b, sizeof b, "%.4f", v);
  std::string s(b);
  s.erase(s.find_last_not_of('0') + 1);
  if (!s.empty() && s.back() == '.') s.pop_back();
  return s;
}

inline std::string median_string(std::vector<double> &vals, const char *na = "NA") {
  return format_median(median_of(vals), na);
}

#endif
