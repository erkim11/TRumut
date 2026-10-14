/*
calling.h: mutation calling of a test sample against references
*/

#ifndef CALLING_H
#define CALLING_H

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

#include "catalog.h"

static const size_t MIN_UNIQUE_READS_PER_INPUT = 10;

// --input-type, shared by cancer/trio/use-pon and matching trcall allele
//   0 = not given, 1 = reads, -1 = assembly
template <class T>
inline bool parse_int_arg(const char* flag, const char* s, T& out) {
  char* end = nullptr;
  errno = 0;
  long long v = std::strtoll(s, &end, 10);
  if (!s[0] || !end || *end || errno) {
    std::cerr << "[ERROR] " << flag << " must be an integer (got '" << s << "')\n";
    return false;
  }
  out = (T)v;
  return true;
}
inline bool parse_double_arg(const char* flag, const char* s, double& out) {
  char* end = nullptr;
  errno = 0;
  double v = std::strtod(s, &end);
  if (!s[0] || !end || *end || errno) {
    std::cerr << "[ERROR] " << flag << " must be a number (got '" << s << "')\n";
    return false;
  }
  out = v;
  return true;
}
inline int parse_input_type(const char* arg) {
  if (std::strcmp(arg, "reads") == 0)    return 1;
  if (std::strcmp(arg, "assembly") == 0) return -1;
  return 0;
}
static const char* INPUT_TYPE_REQUIRED_MSG =
  "--input-type is required: 'reads' (allele output of sequencing reads, aligned or\n"
  "not) or 'assembly' (allele output of haplotype contigs)\n";

static const size_t MIN_READS_INVALID = (size_t)-1;
inline size_t resolve_min_reads(bool assembly, int ploidy, int fixed_min_reads,
                                const char* err_prefix) {
  if (assembly && ploidy == 0) {
    std::cerr << err_prefix << "--ploidy is required with --input-type assembly (1 haploid, 2 diploid, ...)\n";
    return MIN_READS_INVALID;
  }
  if (ploidy > 0) {
    if (!assembly) {
      std::cerr << err_prefix << "--ploidy applies to --input-type assembly only\n";
      return MIN_READS_INVALID;
    }
    if (fixed_min_reads >= 0)
      std::cerr << "[WARNING] --min-reads-per-sample " << fixed_min_reads << " ignored: --ploidy "
                << ploidy << " sets the read floor\n";
    return (size_t)ploidy;
  }
  return fixed_min_reads >= 0 ? (size_t)fixed_min_reads : MIN_UNIQUE_READS_PER_INPUT;
}

inline void die(const std::string& msg) {
  std::cerr << "[ERROR] " << msg << "\n";
  std::exit(1);
}

inline void split_tsv_views(const std::string& line, std::vector<std::string_view>& out) {
  out.clear();
  const char* s = line.data();
  size_t n = line.size();
  size_t i = 0, start = 0;
  while (i <= n) {
    if (i == n || line[i] == '\t') {
      out.emplace_back(s + start, i - start);
      start = i + 1;
    }
    ++i;
  }
}

// parse a field as an integer.
inline int64_t sv_to_i64(std::string_view sv) {
  int64_t x = 0;
  bool neg = false;
  size_t i = 0;
  if (!sv.empty() && sv[0] == '-') { neg = true; i = 1; }
  for (; i < sv.size(); ++i) x = x * 10 + (sv[i] - '0');
  return neg ? -x : x;
}

inline std::string sv_to_string(std::string_view sv) {
  return std::string(sv.begin(), sv.end());
}

// column indices (0-based)
struct ColIdx {
  int chr=-1, region_start=-1, region_end=-1, loci=-1, kmerA=-1, kmerB=-1, allele_lengths=-1;
  int n_motif_kmers=-1;
  int allele_seqs=-1;
  // optional haplotag list.
  int hp=-1;
};

// split a comma-separated list column into items.
inline void split_list(std::string_view field, std::vector<std::string_view>& out) {
  out.clear();
  if (field.empty()) return;
  size_t p = 0;
  while (true) {
    size_t c = field.find(',', p);
    if (c == std::string_view::npos) { out.push_back(field.substr(p)); return; }
    out.push_back(field.substr(p, c - p));
    p = c + 1;
  }
}

struct RowReads {
  size_t n = 0;
  std::vector<std::string_view> lens, motifs, hp, seqs;
  static std::string_view at(const std::vector<std::string_view>& v, size_t i) {
    return i < v.size() ? v[i] : std::string_view(".");
  }
};
inline void expand_row_reads(const std::vector<std::string_view>& f, const ColIdx& ci, RowReads& r) {
  split_list(f[ci.allele_lengths], r.lens);
  r.n = r.lens.size();
  split_list(f[ci.n_motif_kmers], r.motifs);
  if (ci.hp >= 0 && (int)f.size() > ci.hp) split_list(f[ci.hp], r.hp); else r.hp.clear();
  if (ci.allele_seqs >= 0 && (int)f.size() > ci.allele_seqs) split_list(f[ci.allele_seqs], r.seqs); else r.seqs.clear();
}

// parse TSV header and fill ColIdx with 0-based column indices for required fields.
inline bool parse_header_indices(const std::string& header, ColIdx& ci) {
  std::vector<std::string_view> f;
  f.reserve(64);
  split_tsv_views(header, f);
  for (int i = 0; i < (int)f.size(); ++i) {
    const auto name = f[i];
    if (name == "chr") ci.chr = i;
    else if (name == "region_start") ci.region_start = i;
    else if (name == "region_end") ci.region_end = i;
    else if (name == "loci") ci.loci = i;
    else if (name == "kmerA_id") ci.kmerA = i;
    else if (name == "kmerB_id") ci.kmerB = i;
    else if (name == "allele_lengths") ci.allele_lengths = i;
    else if (name == "n_motif_kmers") ci.n_motif_kmers = i;
    else if (name == "allele_seqs") ci.allele_seqs = i;
    else if (name == "HP") ci.hp = i;
  }
  return (ci.chr>=0 && ci.region_start>=0 && ci.region_end>=0 && ci.loci>=0 &&
          ci.kmerA>=0 && ci.kmerB>=0 && ci.allele_lengths>=0 && ci.n_motif_kmers>=0);
}

// create output directory.
inline void ensure_dir_exists(const std::string& dir) {
  if (dir.empty()) die("output_dir is empty");
  std::string cmd = "bash -lc 'mkdir -p " + dir + "'";
  int rc = std::system(cmd.c_str());
  if (rc != 0) die("failed to create output dir (rc=" + std::to_string(rc) + ")");
}

// kmer canonical key -> count
using KmerCntMap = std::unordered_map<uint64_t, uint32_t>;
// kmer canonical key -> sorted list of 0-based start positions
using KmerPosMap = std::unordered_map<uint64_t, std::vector<int32_t>>;

// rolling hash fill: count only.
inline void kmer_fill_counts(KmerCntMap& h, int k, int len, const char* seq) {
  const uint64_t mask = (1ULL << (k * 2)) - 1;
  const int shift = (k - 1) * 2;
  uint64_t x0 = 0, x1 = 0;
  int l = 0;
  for (int i = 0; i < len; ++i) {
    int c = (int)kom_nt4_table[(uint8_t)seq[i]];
    if (c < 4) {
      x0 = (x0 << 2 | (uint64_t)c) & mask;
      x1 = (x1 >> 2) | ((uint64_t)(3 - c) << shift);
      if (++l >= k) ++h[x0 < x1 ? x0 : x1];
    } else { l = 0; x0 = x1 = 0; }
  }
}

// rolling hash fill: record 0-based start positions.
inline void kmer_fill_positions(KmerPosMap& h, int k, int len, const char* seq) {
  const uint64_t mask = (1ULL << (k * 2)) - 1;
  const int shift = (k - 1) * 2;
  uint64_t x0 = 0, x1 = 0;
  int l = 0;
  for (int i = 0; i < len; ++i) {
    int c = (int)kom_nt4_table[(uint8_t)seq[i]];
    if (c < 4) {
      x0 = (x0 << 2 | (uint64_t)c) & mask;
      x1 = (x1 >> 2) | ((uint64_t)(3 - c) << shift);
      if (++l >= k) h[x0 < x1 ? x0 : x1].push_back(i - k + 1);
    } else { l = 0; x0 = x1 = 0; }
  }
}

inline void fill_kmer_counts_for_seq(KmerCntMap& h, const std::string& seq, int k) {
  if ((int)seq.size() < k) return;
  kmer_fill_counts(h, k, (int)seq.size(), seq.data());
}
inline void fill_kmer_pos_map_for_seq(KmerPosMap& h, const std::string& seq, int k) {
  if ((int)seq.size() < k) return;
  kmer_fill_positions(h, k, (int)seq.size(), seq.data());
}

// for each gap between consecutive motif-kmer hit positions in the test sequence:
//   - gap_total[i] = number of non-motif kmer hits in gap i
//   - gap_diff[i]  = how many of those are "excess" in test vs reference (test_count > ref_count)
inline int64_t compute_nonmotif_kmers_in(
    const KmerPosMap& h_test,
    const KmerCntMap& h_ref,
    const std::unordered_set<uint64_t>& motif_kmers,
    int k, int32_t min_gap_kmers, double min_gap_frac)
{
  // 1. diff_keys: kmers where test count exceeds reference count
  std::unordered_set<uint64_t> diff_keys;
  diff_keys.reserve(h_test.size());
  for (const auto& kv : h_test) {
    uint64_t key = kv.first;
    uint32_t t_cnt = (uint32_t)kv.second.size();
    auto it = h_ref.find(key);
    uint32_t r_cnt = it != h_ref.end() ? it->second : 0;
    if (t_cnt > r_cnt) diff_keys.insert(key);
  }

  // 2. collect all motif-kmer hit positions in the test sequence
  std::vector<int32_t> motif_positions;
  for (uint64_t mk : motif_kmers) {
    auto it = h_test.find(mk);
    if (it != h_test.end()) {
      for (int32_t p : it->second) motif_positions.push_back(p);
    }
  }
  std::sort(motif_positions.begin(), motif_positions.end());

  int n_gaps = (int)motif_positions.size() - 1;
  if (n_gaps <= 0) return 0;

  // 3. assign non-motif kmer positions to gaps
  std::vector<int32_t> gap_total(n_gaps, 0);
  std::vector<int32_t> gap_diff(n_gaps, 0);

  for (const auto& kv : h_test) {
    uint64_t key = kv.first;
    if (motif_kmers.count(key)) continue;
    bool in_diff = diff_keys.count(key) > 0;
    for (int32_t pos : kv.second) {
      auto it = std::upper_bound(motif_positions.begin(), motif_positions.end(), pos - k);
      int i = (int)(it - motif_positions.begin()) - 1;
      if (i < 0 || i >= n_gaps) continue;
      int32_t L = motif_positions[i] + k;
      int32_t R = motif_positions[i + 1] - 1;
      if (pos >= L && pos <= R) {
        ++gap_total[i];
        if (in_diff) ++gap_diff[i];
      }
    }
  }

  // 4. sum qualifying gaps
  int64_t result = 0;
  for (int i = 0; i < n_gaps; ++i) {
    int32_t L = motif_positions[i] + k;
    int32_t R = motif_positions[i + 1] - 1;
    if (L > R) continue;
    if (gap_total[i] < min_gap_kmers) continue;
    if ((double)gap_diff[i] / (double)gap_total[i] >= min_gap_frac) result += gap_diff[i];
  }
  return result;
}

struct HpCounts {
  int32_t hp1 = 0;
  int32_t hp2 = 0;
  int32_t na  = 0;
};

struct GroupData {
  // metadata from test group
  std::string chr;
  int64_t min_start = std::numeric_limits<int64_t>::max();
  int64_t max_end   = std::numeric_limits<int64_t>::min();
  std::string kmerA_id;
  std::string kmerB_id;
  std::string loci;   // "start-end:motif" tokens

  // frequencies
  std::unordered_map<int32_t,int32_t> test_freq; // test length -> count
  std::unordered_map<int32_t,int32_t> ref_freq;  // reference length -> count
  std::unordered_map<int32_t,HpCounts> test_hp_freq;
  size_t n_test_reads = 0;
  std::vector<size_t> n_ref_reads;

  std::unordered_map<int32_t, std::vector<double>> test_motif_by_len;
  std::unordered_map<int32_t, std::vector<double>> ref_motif_by_len;

  // every reference read's length
  std::vector<int32_t> ref_distances;
  std::vector<std::string> test_seqs;
  std::vector<std::string> ref_seqs;
  std::vector<int32_t> test_seq_lens;   // parallel to test_seqs
  std::vector<int32_t> ref_seq_lens;    // parallel to ref_seqs
};

// update group metadata and test length frequencies from one test row
inline void update_test_meta_and_freq(GroupData& g, const std::vector<std::string_view>& f, const ColIdx& ci, RowReads& rr) {
  if (g.chr.empty())  g.chr  = sv_to_string(f[ci.chr]);
  if (g.loci.empty()) g.loci = sv_to_string(f[ci.loci]);
  int64_t st = sv_to_i64(f[ci.region_start]);
  int64_t en = sv_to_i64(f[ci.region_end]);
  if (st < g.min_start) g.min_start = st;
  if (en > g.max_end)   g.max_end   = en;

  expand_row_reads(f, ci, rr);
  for (size_t i = 0; i < rr.n; ++i) {
    const int32_t L = (int32_t)sv_to_i64(rr.lens[i]);
    ++g.test_freq[L];
    ++g.n_test_reads;
    g.test_motif_by_len[L].push_back(parse_motif_count(RowReads::at(rr.motifs, i)));
    const std::string_view h = RowReads::at(rr.hp, i);
    if      (h == "1") ++g.test_hp_freq[L].hp1;
    else if (h == "2") ++g.test_hp_freq[L].hp2;
    else               ++g.test_hp_freq[L].na;
    if (ci.allele_seqs >= 0) {
      std::string_view sq = RowReads::at(rr.seqs, i);
      g.test_seqs.emplace_back(sq == "." ? std::string_view() : sq);
      g.test_seq_lens.push_back(L);
    }
  }
}

// update group reference length frequencies from one reference row
inline void update_ref_freq(GroupData& g, const std::vector<std::string_view>& f, const ColIdx& ci, RowReads& rr, size_t file_idx) {
  expand_row_reads(f, ci, rr);
  for (size_t i = 0; i < rr.n; ++i) {
    const int32_t L = (int32_t)sv_to_i64(rr.lens[i]);
    ++g.ref_freq[L];
    g.ref_distances.push_back(L);
    g.ref_motif_by_len[L].push_back(parse_motif_count(RowReads::at(rr.motifs, i)));
    if (ci.allele_seqs >= 0) {
      std::string_view sq = RowReads::at(rr.seqs, i);
      g.ref_seqs.emplace_back(sq == "." ? std::string_view() : sq);
      g.ref_seq_lens.push_back(L);
    }
    ++g.n_ref_reads[file_idx];
  }
}

// return minimum absolute difference from L to any element in S.
inline int32_t nearest_diff_to_set(int32_t L, const std::vector<int32_t>& S) {
  int32_t best = std::numeric_limits<int32_t>::max();
  for (int32_t x : S) {
    int32_t d = std::abs(L - x);
    if (d < best) best = d;
  }
  return best;
}

inline int32_t signed_nearest_diff_to_set(int32_t L, const std::vector<int32_t>& S) {
  int32_t best = std::numeric_limits<int32_t>::max();
  for (int32_t x : S) {
    int32_t d = L - x;
    if (std::abs(d) < std::abs(best) || (std::abs(d) == std::abs(best) && d > best)) best = d;
  }
  return best;
}

// a candidate deviating length with its supporting-read count.
struct LenCand { int32_t len; int32_t cnt; };

inline std::vector<int32_t> greedy_cluster_lengths(
    std::vector<LenCand> cands,
    int32_t threshold)
{
  if (cands.empty()) return {};

  auto pick_best = [](const std::vector<LenCand>& v)->int32_t{
    int32_t bestL = v[0].len, bestC = v[0].cnt;
    for (const auto& c : v) {
      if (c.cnt > bestC || (c.cnt == bestC && c.len > bestL)) {
        bestC = c.cnt; bestL = c.len;
      }
    }
    return bestL;
  };

  std::vector<int32_t> deviating;
  deviating.reserve(4);

  std::vector<LenCand> remaining = std::move(cands);
  while (!remaining.empty()) {
    int32_t chosen = pick_best(remaining);
    deviating.push_back(chosen);

    std::vector<LenCand> next;
    next.reserve(remaining.size());
    for (const auto& c : remaining) {
      if (c.len == chosen) continue;
      if (nearest_diff_to_set(c.len, deviating) > threshold) next.push_back(c);
    }
    remaining.swap(next);
  }
  return deviating;
}

// return test lengths that deviate from every reference length by at least min_offset bp.
inline std::vector<int32_t> find_deviating_lengths(
    const std::unordered_map<int32_t,int32_t>& test_freq,
    const std::vector<int32_t>& ref_lens,
    int32_t threshold,
    int32_t min_offset)
{
  if (ref_lens.empty()) return {};
  std::vector<LenCand> cands;
  cands.reserve(test_freq.size());

  for (const auto& kv : test_freq) {
    int32_t L = kv.first, cnt = kv.second;
    if (nearest_diff_to_set(L, ref_lens) >= min_offset) cands.push_back({L, cnt});
  }
  return greedy_cluster_lengths(std::move(cands), threshold);
}

inline std::string parent_dir(const std::string& path) {
  size_t slash = path.find_last_of('/');
  if (slash == std::string::npos) return ".";
  if (slash == 0) return "/";
  return path.substr(0, slash);
}

inline std::string base_name(const std::string& path) {
  size_t slash = path.find_last_of('/');
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

struct CallLabels {
  const char* test_label;
  const char* ref_label;
};

struct RowFilters {
  int32_t min_supporting_reads = 7;
  double  min_motif_diff = 3.0;
  int32_t min_subclone_reads = 4;
  int32_t min_gap_kmers = 10;
  double  min_gap_frac = 0.95;
  bool keep(bool assembly, int32_t n_reads_in_thres, double motif_test, double motif_ref) const {
    if (!assembly && min_supporting_reads >= 0 && n_reads_in_thres < min_supporting_reads) return false;
    if (min_motif_diff >= 0 && !std::isnan(motif_test) && !std::isnan(motif_ref) &&
        std::fabs(motif_test - motif_ref) <= min_motif_diff) return false;
    return true;
  }
};

inline std::string alleles_header(const CallLabels& labels, bool with_phase = false,
                                  bool assembly = false) {
  (void)labels;
  std::string h = "chr\tregion_start\tregion_end\tkmerA_id\tkmerB_id\tloci\tallele_length\toffset"
                  "\tn_motif_kmers_ref_median\tn_motif_kmers_test_median\tn_nonmotif_kmers_in";
  if (!assembly) h += "\tn_reads\tn_reads_in_thres";
  h += "\tevent_type";
  if (with_phase) h += "\tphase_counts";
  h += "\n";
  return h;
}

// number of test reads in the group whose length is exactly `len`.
inline int32_t exact_read_count(const std::unordered_map<int32_t,int32_t>& test_freq,
                                int32_t len) {
  auto it = test_freq.find(len);
  return it == test_freq.end() ? 0 : it->second;
}

// number of test reads in the group whose length is within `threshold` of `len` (inclusive)
inline int32_t nearby_read_count(const std::unordered_map<int32_t,int32_t>& test_freq,
                                 int32_t len, int32_t threshold) {
  int32_t n = 0;
  for (const auto& kv : test_freq) {
    if (std::abs(kv.first - len) <= threshold) n += kv.second;
  }
  return n;
}

// dominant haplotype phase of the test reads within `threshold` of `len`
enum class DomPhase { HP1, HP2, Tie, None };

inline DomPhase dominant_phase_for_length(
    const std::unordered_map<int32_t,HpCounts>& test_hp_freq,
    int32_t len, int32_t threshold) {
  int32_t n1 = 0, n2 = 0;
  for (const auto& kv : test_hp_freq) {
    if (std::abs(kv.first - len) <= threshold) {
      n1 += kv.second.hp1;
      n2 += kv.second.hp2;
    }
  }
  if (n1 == 0 && n2 == 0) return DomPhase::None;
  if (n1 == n2)           return DomPhase::Tie;
  return n1 > n2 ? DomPhase::HP1 : DomPhase::HP2;
}

inline std::unordered_map<int32_t,int32_t> pooled_exact_freq(
    const std::unordered_map<int32_t,HpCounts>& test_hp_freq) {
  std::unordered_map<int32_t,int32_t> out;
  out.reserve(test_hp_freq.size());
  for (const auto& kv : test_hp_freq) {
    const int32_t n = kv.second.hp1 + kv.second.hp2 + kv.second.na;
    if (n > 0) out[kv.first] = n;
  }
  return out;
}

inline bool distant_same_phase_cluster(
    const std::unordered_map<int32_t,int32_t>& freq,
    int32_t len, int32_t threshold, int32_t min_reads) {
  std::vector<std::pair<int32_t,int32_t>> pts;   // (length, phase read count)
  for (const auto& kv : freq)
    if (kv.second > 0) pts.emplace_back(kv.first, kv.second);
  pts.emplace_back(len, 0);
  std::sort(pts.begin(), pts.end());
  int32_t cluster_reads = 0;
  bool cluster_has_len = false;
  for (size_t i = 0; i < pts.size(); ++i) {
    if (i > 0 && pts[i].first - pts[i - 1].first > threshold) {
      if (!cluster_has_len && cluster_reads >= min_reads) return true;
      cluster_reads = 0;
      cluster_has_len = false;
    }
    cluster_reads += pts[i].second;
    if (pts[i].first == len) cluster_has_len = true;
  }
  return !cluster_has_len && cluster_reads >= min_reads;
}

inline bool is_private_length(const std::unordered_map<int32_t,int32_t>& freq,
                              int32_t len, int32_t threshold) {
  auto it = freq.find(len);
  const int32_t f = (it == freq.end()) ? 0 : it->second;
  for (const auto& kv : freq) {
    if (std::abs(kv.first - len) <= threshold && kv.second > f) return true;
  }
  return false;
}

inline std::unordered_map<int32_t,int32_t> phase_exact_freq(
    const std::unordered_map<int32_t,HpCounts>& test_hp_freq, DomPhase hp) {
  std::unordered_map<int32_t,int32_t> out;
  out.reserve(test_hp_freq.size());
  for (const auto& kv : test_hp_freq) {
    const int32_t n = (hp == DomPhase::HP1) ? kv.second.hp1 : kv.second.hp2;
    if (n > 0) out[kv.first] = n;
  }
  return out;
}

inline int32_t other_phase_reads_in_window(
    const std::unordered_map<int32_t,HpCounts>& test_hp_freq,
    int32_t len, int32_t threshold, DomPhase hp) {
  int32_t n = 0;
  for (const auto& kv : test_hp_freq)
    if (std::abs(kv.first - len) <= threshold)
      n += (hp == DomPhase::HP1) ? kv.second.hp2 : kv.second.hp1;
  return n;
}

inline void add_untagged_in_window(
    std::unordered_map<int32_t,int32_t>& freq,
    const std::unordered_map<int32_t,HpCounts>& test_hp_freq,
    int32_t len, int32_t threshold) {
  for (const auto& kv : test_hp_freq)
    if (kv.second.na > 0 && std::abs(kv.first - len) <= threshold)
      freq[kv.first] += kv.second.na;
}

using LengthPredicate = std::function<bool(int32_t)>;

inline bool window_mode(const std::unordered_map<int32_t,int32_t>& freq,
                        int32_t len, int32_t threshold, int32_t& mode,
                        const LengthPredicate& eligible = {}) {
  bool found = false;
  int32_t best_cnt = 0;
  for (const auto& kv : freq) {
    if (kv.second <= 0 || std::abs(kv.first - len) > threshold) continue;
    if (eligible && !eligible(kv.first)) continue;
    const bool better =
        !found || kv.second > best_cnt ||
        (kv.second == best_cnt &&
         (std::abs(kv.first - len) < std::abs(mode - len) ||
          (std::abs(kv.first - len) == std::abs(mode - len) && kv.first > mode)));
    if (better) { found = true; best_cnt = kv.second; mode = kv.first; }
  }
  return found;
}

inline void warn_unphased_once() {
  static std::atomic<bool> warned{false};
  if (warned.exchange(true)) return;
  std::cerr << "[WARNING] some loci have no HP-tagged read in the window; their event_type is phase-free (private_unphased / not_private_unphased)\n";
}

inline const char* event_type_for_length(
    const std::unordered_map<int32_t,HpCounts>& test_hp_freq,
    int32_t len, int32_t threshold,
    const LengthPredicate& is_deviating = {},
    int32_t min_subclone_reads = 4)
{
  const DomPhase cand = dominant_phase_for_length(test_hp_freq, len, threshold);
  {
    auto it = test_hp_freq.find(len);
    if (it != test_hp_freq.end() &&
        it->second.hp1 + it->second.hp2 + it->second.na == 1)
      return cand == DomPhase::None ? "private_unphased" : "private_phased";
  }
  if (cand == DomPhase::None) {
    if (test_hp_freq.empty()) return "N/A";
    warn_unphased_once();
    std::unordered_map<int32_t,int32_t> freq = pooled_exact_freq(test_hp_freq);
    int32_t mode = len;
    if (!window_mode(freq, len, threshold, mode, is_deviating)) return "private_unphased";
    if (is_private_length(freq, mode, threshold)) return "private_unphased";
    return "not_private_unphased";
  }
  bool any_peak = false;
  const DomPhase phases[2] = { DomPhase::HP1, DomPhase::HP2 };
  for (DomPhase hp : phases) {
    if (cand != DomPhase::Tie && cand != hp) continue;
    std::unordered_map<int32_t,int32_t> freq = phase_exact_freq(test_hp_freq, hp);
    if (other_phase_reads_in_window(test_hp_freq, len, threshold, hp) == 0)
      add_untagged_in_window(freq, test_hp_freq, len, threshold);
    int32_t mode = len;
    if (!window_mode(freq, len, threshold, mode, is_deviating)) continue;
    if (is_private_length(freq, mode, threshold)) continue;
    any_peak = true;
    if (distant_same_phase_cluster(freq, mode, threshold, min_subclone_reads))
      return "subclonal";
  }
  return any_peak ? "clonal" : "private_phased";
}

// haplotag split of the same reads nearby_read_count counts: of the test reads within `threshold` of `len`.
inline std::string phase_counts_for_length(
    const std::unordered_map<int32_t,HpCounts>& test_hp_freq,
    int32_t len, int32_t threshold) {
  int32_t n1 = 0, n2 = 0, nna = 0;
  for (const auto& kv : test_hp_freq) {
    if (std::abs(kv.first - len) <= threshold) {
      n1 += kv.second.hp1; n2 += kv.second.hp2; nna += kv.second.na;
    }
  }
  return std::to_string(n1) + ":" + std::to_string(n2) + ":" + std::to_string(nna);
}

// median n_motif_kmers of the reads at length L.
inline double motif_median_value(const std::unordered_map<int32_t, std::vector<double>>& by_len,
                                 int32_t L) {
  auto it = by_len.find(L);
  if (it == by_len.end()) return std::numeric_limits<double>::quiet_NaN();
  std::vector<double> vals = it->second;
  return median_of(vals);
}

inline int32_t nearest_ref_length(int32_t d, const std::unordered_map<int32_t,int32_t>& ref_freq) {
  int32_t best = 0;
  int64_t best_diff = std::numeric_limits<int64_t>::max();
  for (const auto& kv : ref_freq) {
    int64_t diff = std::llabs((int64_t)kv.first - d);
    if (diff < best_diff || (diff == best_diff && kv.first < best)) { best_diff = diff; best = kv.first; }
  }
  return best;
}

inline void write_alleles_for_group(
    std::ofstream& out,
    const GroupData& g,
    const std::vector<int32_t>& deviating,
    const ColIdx& read_ci,
    const std::unordered_set<uint64_t>& union_motif_kmers,
    const CallLabels& labels,
    int32_t n_reads_thres,
    int32_t min_offset,
    const RowFilters& filters,
    bool assembly = false)
{
  const std::string& chr = g.chr;
  (void)labels;

  const int K = MOTIF_KMER_K;
  const bool have_seqs = (read_ci.allele_seqs >= 0) && !g.test_seqs.empty() && !g.ref_seqs.empty();
  const bool with_phase = (read_ci.hp >= 0);

  {
    for (int32_t d : deviating) {
      int32_t offset_val = signed_nearest_diff_to_set(d, g.ref_distances);
      const char* event_type = assembly
          ? "clonal"
          : event_type_for_length(g.test_hp_freq, d, n_reads_thres,
                [&](int32_t L) { return nearest_diff_to_set(L, g.ref_distances) >= min_offset; },
                filters.min_subclone_reads);

      const int32_t num_exact_val  = exact_read_count(g.test_freq, d);
      const int32_t num_nearby_val = nearby_read_count(g.test_freq, d, n_reads_thres);
      const std::string phase = with_phase
          ? phase_counts_for_length(g.test_hp_freq, d, n_reads_thres) : std::string();

      const double motif_test_v = motif_median_value(g.test_motif_by_len, d);
      const double motif_ref_v  = g.ref_freq.empty() ? std::numeric_limits<double>::quiet_NaN()
          : motif_median_value(g.ref_motif_by_len, nearest_ref_length(d, g.ref_freq));

      if (!filters.keep(assembly, num_nearby_val, motif_test_v, motif_ref_v)) continue;

      const std::string motif_test = format_median(motif_test_v, "nan");
      const std::string motif_ref  = format_median(motif_ref_v, "nan");

      std::string nonmotif_in = "nan";
      if (have_seqs && !union_motif_kmers.empty()) {
        const std::string* test_seq = nullptr;
        for (size_t i = 0; i < g.test_seqs.size(); ++i)
          if (g.test_seq_lens[i] == d) { test_seq = &g.test_seqs[i]; break; }
        const std::string* ref_seq = nullptr;
        int32_t best_diff = std::numeric_limits<int32_t>::max();
        for (size_t i = 0; i < g.ref_seqs.size(); ++i) {
          int32_t diff = std::abs(g.ref_seq_lens[i] - d);
          if (diff < best_diff) { best_diff = diff; ref_seq = &g.ref_seqs[i]; }
        }
        if (test_seq && ref_seq) {
          KmerCntMap h_ref;
          KmerPosMap h_test;
          fill_kmer_counts_for_seq(h_ref, *ref_seq, K);
          fill_kmer_pos_map_for_seq(h_test, *test_seq, K);
          nonmotif_in = std::to_string(compute_nonmotif_kmers_in(h_test, h_ref, union_motif_kmers, K, filters.min_gap_kmers, filters.min_gap_frac));
        }
      }

      out << chr << '\t' << g.min_start << '\t' << g.max_end << '\t'
          << g.kmerA_id << '\t' << g.kmerB_id << '\t' << g.loci
          << '\t' << d << '\t' << offset_val
          << '\t' << motif_ref << '\t' << motif_test
          << '\t' << nonmotif_in;
      if (!assembly) out << '\t' << num_exact_val << '\t' << num_nearby_val;
      out << '\t' << event_type;
      if (with_phase) out << '\t' << phase;
      out << '\n';
    }
  }
}

inline int cmp_key(const std::string& aA, const std::string& aB,
                   const std::string& bA, const std::string& bB) {
  if (aA < bA) return -1;
  if (aA > bA) return  1;
  if (aB < bB) return -1;
  if (aB > bB) return  1;
  return 0;
}

struct ChrRange {
  std::string chr;
  int64_t begin = 0;
  int64_t end   = 0;
};

struct InputIndex {
  std::string path;
  bool is_temp = false;
  std::string header;
  ColIdx ci;
  std::vector<ChrRange> ranges;   // in file order
  const ChrRange* find(const std::string& chr) const {
    for (const auto& r : ranges) if (r.chr == chr) return &r;
    return nullptr;
  }
};

inline bool is_gzip_file(const std::string& path) {
  std::FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) return false;
  unsigned char m[2];
  size_t n = std::fread(m, 1, 2, f);
  std::fclose(f);
  return n == 2 && m[0] == 0x1f && m[1] == 0x8b;
}

struct LineReader {
  gzFile g = nullptr;
  bool open(const std::string& path) {
    g = gzopen(path.c_str(), "rb");
    return g != nullptr;
  }
  explicit operator bool() const { return g != nullptr; }
  bool getline(std::string& line) {
    line.clear();
    if (!g) return false;
    char buf[1 << 16];
    for (;;) {
      if (!gzgets(g, buf, sizeof buf)) return !line.empty();
      size_t n = std::strlen(buf);
      if (n && buf[n - 1] == '\n') { line.append(buf, n - 1); return true; }
      line.append(buf, n);
    }
  }
  ~LineReader() { if (g) gzclose(g); }
  LineReader() = default;
  LineReader(const LineReader&) = delete;
  LineReader& operator=(const LineReader&) = delete;
};

// decompress a gzipped input into outdir as a plain file.
inline std::string gunzip_to_plain(const std::string& in_path, const std::string& tag,
                                   const std::string& outdir) {
  std::string base = in_path;
  {
    size_t slash = base.find_last_of('/');
    if (slash != std::string::npos) base = base.substr(slash + 1);
  }
  std::string out_path = outdir;
  if (!out_path.empty() && out_path.back() != '/') out_path.push_back('/');
  out_path += base + "." + tag + ".plain.tsv";

  gzFile g = gzopen(in_path.c_str(), "rb");
  if (!g) die("cannot open input: " + in_path);
  std::FILE* f = std::fopen(out_path.c_str(), "wb");
  if (!f) { gzclose(g); die("cannot write decompressed copy: " + out_path); }
  std::vector<char> buf(1 << 20);
  int n;
  while ((n = gzread(g, buf.data(), (unsigned)buf.size())) > 0)
    if (std::fwrite(buf.data(), 1, (size_t)n, f) != (size_t)n) die("write failed: " + out_path);
  int err = 0;
  const char* msg = gzerror(g, &err);
  if (n < 0 || err == Z_ERRNO) die("decompressing " + in_path + " failed: " + (msg ? msg : ""));
  gzclose(g);
  if (std::fclose(f) != 0) die("write failed: " + out_path);
  return out_path;
}

inline bool scan_and_index(const std::string& path, InputIndex& ix) {
  std::ifstream in(path);
  if (!in) die("cannot open input: " + path);
  ix.path = path;
  ix.ranges.clear();
  if (!std::getline(in, ix.header)) die("empty input: " + path);
  if (!parse_header_indices(ix.header, ix.ci)) die("missing required columns in header: " + path);

  int64_t off = (int64_t)ix.header.size() + 1;
  std::string line;
  std::vector<std::string_view> f;
  f.reserve(64);
  std::set<std::string> closed;   // chromosomes whose block already ended
  std::string prevA, prevB;
  bool have_prev = false;
  const int need = std::max({ix.ci.chr, ix.ci.kmerA, ix.ci.kmerB});

  while (std::getline(in, line)) {
    const int64_t next = off + (int64_t)line.size() + 1;
    if (!line.empty()) {
      split_tsv_views(line, f);
      if ((int)f.size() > need) {
        std::string chr = sv_to_string(f[ix.ci.chr]);
        if (ix.ranges.empty() || ix.ranges.back().chr != chr) {
          if (!ix.ranges.empty()) closed.insert(ix.ranges.back().chr);
          if (closed.count(chr)) return false;   // chromosome split into two blocks
          ix.ranges.push_back({chr, off, next});
          have_prev = false;
        } else {
          ix.ranges.back().end = next;
        }
        std::string a = sv_to_string(f[ix.ci.kmerA]);
        std::string b = sv_to_string(f[ix.ci.kmerB]);
        if (have_prev && cmp_key(prevA, prevB, a, b) > 0) return false;
        prevA = std::move(a); prevB = std::move(b);
        have_prev = true;
      }
    }
    off = next;
  }
  return true;
}

inline std::string parent_sort_whole_file(const std::string& in_path, const ColIdx& ci,
                                          const std::string& tag, const std::string& outdir) {
  std::string base = in_path;
  {
    size_t slash = base.find_last_of('/');
    if (slash != std::string::npos) base = base.substr(slash + 1);
  }
  std::string out_path = outdir;
  if (!out_path.empty() && out_path.back() != '/') out_path.push_back('/');
  out_path += base + "." + tag + ".sorted.tsv";

  auto key = [](int c) { std::string n = std::to_string(c + 1); return " -k" + n + "," + n; };
  std::string cmd;
  cmd += "bash -lc 'set -euo pipefail; ";
  cmd += "head -n 1 " + in_path + " > " + out_path + "; ";
  cmd += "tail -n +2 " + in_path + " | ";
  cmd += "LC_ALL=C sort -t\"";
  cmd += '\t';
  cmd += "\"" + key(ci.chr) + key(ci.kmerA) + key(ci.kmerB);
  cmd += " -S 1G -T " + outdir;
  cmd += " >> " + out_path + ";'";
  int rc = std::system(cmd.c_str());
  if (rc != 0) die("external sort failed (rc=" + std::to_string(rc) + "). Command was: " + cmd);

  return out_path;
}

inline void build_input_index(const std::string& path, const std::string& tag,
                              const std::string& outdir, InputIndex& ix) {
  std::string src = path;
  bool src_is_temp = false;
  if (is_gzip_file(path)) {
    std::fprintf(stderr, "[INFO] %s is gzipped; decompressing a copy\n", path.c_str());
    src = gunzip_to_plain(path, tag, outdir);
    src_is_temp = true;
  }
  if (scan_and_index(src, ix)) { ix.is_temp = src_is_temp; return; }
  std::fprintf(stderr, "[INFO] %s is not sorted by (chr, kmerA_id, kmerB_id); sorting a copy\n", path.c_str());
  const std::string sorted = parent_sort_whole_file(src, ix.ci, tag, outdir);
  if (src_is_temp) std::remove(src.c_str());
  if (!scan_and_index(sorted, ix)) die("sorting did not produce a sorted file: " + sorted);
  ix.is_temp = true;
}

struct SliceReader {
  std::ifstream in;
  int64_t pos = 0, end = 0;
  void open(const InputIndex& ix, const std::string& chr) {
    const ChrRange* r = ix.find(chr);
    if (!r) return;
    in.open(ix.path);
    if (!in) die("cannot open input: " + ix.path);
    in.seekg(r->begin);
    pos = r->begin;
    end = r->end;
  }
  bool next(std::string& line) {
    if (pos >= end) return false;
    if (!std::getline(in, line)) return false;
    pos += (int64_t)line.size() + 1;
    return true;
  }
};

inline int run_deviation_calling_one_chr(
    const InputIndex& test_ix,
    const std::vector<InputIndex>& ref_ixs,
    const std::string& alleles_path,
    const std::string& chr,
    int32_t fixed_clustering_thres,
    int32_t fixed_n_reads_thres,
    int32_t fixed_min_offset,
    bool with_header,
    const CallLabels& labels,
    const RowFilters& filters,
    size_t min_reads,
    bool assembly)
{
  std::ofstream aout(alleles_path);
  if (!aout) die("cannot write: " + alleles_path);
  if (with_header) aout << alleles_header(labels, test_ix.ci.hp >= 0, assembly);

  SliceReader tin;
  tin.open(test_ix, chr);
  const ColIdx& tci = test_ix.ci;

  // one stream per reference file.
  struct RefStream {
    SliceReader in;
    ColIdx ci;
    std::string line;
    bool have = false;
  };
  std::vector<RefStream> refs(ref_ixs.size());
  for (size_t i = 0; i < ref_ixs.size(); ++i) {
    refs[i].in.open(ref_ixs[i], chr);
    refs[i].ci = ref_ixs[i].ci;
    refs[i].have = refs[i].in.next(refs[i].line);
  }

  std::string tline;
  std::vector<std::string_view> tf, rf;
  tf.reserve(64); rf.reserve(64);
  RowReads rr;

  bool thave = tin.next(tline);

  auto read_key = [&](const std::string& line, std::vector<std::string_view>& f, const ColIdx& ci,
                      std::string& kA, std::string& kB)->bool {
    if (line.empty()) return false;
    split_tsv_views(line, f);
    const int max_col = std::max({ci.kmerA, ci.kmerB});
    if ((int)f.size() <= max_col) return false;
    kA.assign(f[ci.kmerA].begin(), f[ci.kmerA].end());
    kB.assign(f[ci.kmerB].begin(), f[ci.kmerB].end());
    return true;
  };

  const int32_t clustering_thres = fixed_clustering_thres >= 0 ? fixed_clustering_thres : 11;
  const int32_t n_reads_thres = fixed_n_reads_thres >= 0 ? fixed_n_reads_thres : 9;
  const int32_t min_offset = fixed_min_offset >= 0 ? fixed_min_offset : clustering_thres + 1;

  while (thave) {
    std::string curA, curB;
    if (!read_key(tline, tf, tci, curA, curB)) {
      thave = tin.next(tline);
      continue;
    }

    GroupData g;
    g.kmerA_id = curA;
    g.kmerB_id = curB;
    g.n_ref_reads.assign(refs.size(), 0);

    for (size_t ri = 0; ri < refs.size(); ++ri) {
      auto& r = refs[ri];
      // advance this reference until ref_key >= (curA, curB).
      while (r.have) {
        std::string a, b;
        if (!read_key(r.line, rf, r.ci, a, b)) {
          r.have = r.in.next(r.line);
          continue;
        }
        if (cmp_key(a, b, curA, curB) < 0) {
          r.have = r.in.next(r.line);
          continue;
        }
        break;
      }

      while (r.have) {
        std::string a, b;
        if (!read_key(r.line, rf, r.ci, a, b)) {
          r.have = r.in.next(r.line);
          continue;
        }
        if (a != curA || b != curB) break;
        update_ref_freq(g, rf, r.ci, rr, ri);
        r.have = r.in.next(r.line);
      }
    }

    while (thave) {
      std::string a, b;
      if (!read_key(tline, tf, tci, a, b)) {
        thave = tin.next(tline);
        continue;
      }
      if (a != curA || b != curB) break;
      update_test_meta_and_freq(g, tf, tci, rr);
      thave = tin.next(tline);
    }

    // --min-reads-per-sample, default 10
    if (g.n_test_reads < min_reads) continue; // too few test reads for this loci

    if (g.ref_freq.empty()) continue; // no reference data to compare against

    bool ref_underpowered = false;
    for (size_t n : g.n_ref_reads) {
      if (n > 0 && n < min_reads) { ref_underpowered = true; break; }
    }
    if (ref_underpowered) continue;

    // compare test lengths against every reference length in the group
    std::vector<int32_t> ref_lens;
    ref_lens.reserve(g.ref_freq.size());
    for (const auto& kv : g.ref_freq) ref_lens.push_back(kv.first);

    std::vector<int32_t> deviating = find_deviating_lengths(g.test_freq, ref_lens, clustering_thres, min_offset);
    if (deviating.empty()) continue;

    if (g.min_start == std::numeric_limits<int64_t>::max() ||
        g.max_end   == std::numeric_limits<int64_t>::min() ||
        g.chr.empty()) continue; // malformed group

    const std::unordered_set<uint64_t> union_motif_kmers =
        (tci.allele_seqs >= 0 && !g.test_seqs.empty() && !g.ref_seqs.empty())
            ? find_loci_motif_kmers(g.loci) : std::unordered_set<uint64_t>{};

    write_alleles_for_group(aout, g, deviating, tci, union_motif_kmers, labels, n_reads_thres, min_offset, filters, assembly);
  }

  return 0;
}

inline int run_deviation_calling(
    const std::string& test_path,
    const std::vector<std::string>& ref_paths,
    const std::string& out_path,
    const std::string& target_chr,
    int32_t fixed_clustering_thres,
    int32_t fixed_n_reads_thres,
    int32_t fixed_min_offset,
    const CallLabels& labels,
    const RowFilters& filters,
    size_t min_reads,
    bool assembly,
    int n_jobs_max = 3)
{
  if (out_path.empty()) die("output path is empty");
  const std::string outdir = parent_dir(out_path);
  ensure_dir_exists(outdir);

  InputIndex test_ix;
  build_input_index(test_path, labels.test_label, outdir, test_ix);
  std::vector<InputIndex> ref_ixs(ref_paths.size());
  for (size_t i = 0; i < ref_paths.size(); ++i) {
    std::string tag = labels.ref_label;
    if (ref_paths.size() > 1) tag += std::to_string(i + 1);
    build_input_index(ref_paths[i], tag, outdir, ref_ixs[i]);
  }
  auto cleanup_temps = [&]() {
    if (test_ix.is_temp) std::remove(test_ix.path.c_str());
    for (const auto& r : ref_ixs) if (r.is_temp) std::remove(r.path.c_str());
  };

  if (test_ix.ranges.empty()) die("no data rows in: " + test_path);
  std::vector<std::string> chrs;
  if (!target_chr.empty()) {
    chrs.push_back(target_chr);
  } else {
    chrs.reserve(test_ix.ranges.size());
    for (const auto& r : test_ix.ranges) chrs.push_back(r.chr);
  }

  if (chrs.size() == 1) {
    int rc = run_deviation_calling_one_chr(test_ix, ref_ixs, out_path, chrs[0],
                                           fixed_clustering_thres, fixed_n_reads_thres, fixed_min_offset,
                                           /*with_header=*/true, labels, filters, min_reads, assembly);
    cleanup_temps();
    return rc;
  }

  const std::string out_base = base_name(out_path);
  std::vector<std::string> pieces(chrs.size());
  for (size_t i = 0; i < chrs.size(); ++i)
    pieces[i] = outdir + "/." + out_base + "." + chrs[i] + ".partial";

  const int n_jobs = (int)std::min(chrs.size(), (size_t)std::max(1, n_jobs_max));
  std::fprintf(stderr, "[INFO] running %zu chromosomes, %d at a time; use -j to change\n", chrs.size(), n_jobs);

  size_t next = 0;
  int running = 0, n_failed = 0, n_done = 0;
  std::vector<pid_t> pids(chrs.size(), -1);
  std::vector<char> ok(chrs.size(), 0);

  while (next < chrs.size() || running > 0) {
    while (running < n_jobs && next < chrs.size()) {
      std::fflush(nullptr);
      pid_t pid = fork();
      if (pid < 0) {
        std::fprintf(stderr, "[ERROR] cannot fork a worker for %s\n", chrs[next].c_str());
        ++n_failed; ++next;
        continue;
      }
      if (pid == 0)
        _exit(run_deviation_calling_one_chr(test_ix, ref_ixs, pieces[next], chrs[next],
                                            fixed_clustering_thres, fixed_n_reads_thres, fixed_min_offset,
                                            /*with_header=*/false, labels, filters, min_reads, assembly));
      pids[next] = pid;
      ++running; ++next;
    }

    int st = 0;
    pid_t done_pid = waitpid(-1, &st, 0);
    if (done_pid < 0) {
      if (errno == EINTR) continue;
      break;
    }
    for (size_t i = 0; i < pids.size(); ++i) {
      if (pids[i] != done_pid) continue;
      pids[i] = -1;
      --running; ++n_done;
      const int rc = WIFEXITED(st) ? WEXITSTATUS(st) : -1;
      if (rc == 0) {
        ok[i] = 1;
        std::fprintf(stderr, "  [%2d/%zu] %s done\n", n_done, chrs.size(), chrs[i].c_str());
      } else {
        ++n_failed;
        std::fprintf(stderr, "  [%2d/%zu] %s FAILED (rc=%d)\n", n_done, chrs.size(), chrs[i].c_str(), rc);
      }
      break;
    }
  }

  cleanup_temps();

  {
    std::ofstream fout(out_path);
    if (!fout) die("cannot write: " + out_path);
    fout << alleles_header(labels, test_ix.ci.hp >= 0, assembly);
    for (size_t i = 0; i < chrs.size(); ++i) {
      if (ok[i]) {
        std::ifstream pin(pieces[i], std::ios::binary);
        if (pin && pin.peek() != std::char_traits<char>::eof()) fout << pin.rdbuf();
      }
      std::remove(pieces[i].c_str());
    }
    if (!fout) die("writing failed: " + out_path);
  }

  if (n_failed)
    std::fprintf(stderr, "[ERROR] %zu of %zu chromosomes failed; %s is incomplete (%zu succeeded), see messages above\n",
                 (size_t)n_failed, chrs.size(), out_path.c_str(), chrs.size() - (size_t)n_failed);
  return n_failed ? 1 : 0;
}

#endif
