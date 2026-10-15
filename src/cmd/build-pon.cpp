/*
build-pon.cpp -- `trumut build-pon`: builds a panel of normals for `trumut allele` from your own tandem repeat catalog
*/

#include <fcntl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "lib/calling.h"

namespace fs = std::filesystem;

namespace {

inline bool ends_with(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

// fast buffered line reader over one byte range [start, end) of a file;
// end == -1 reads to EOF
class FastLineReader {
public:
  FastLineReader(const std::string& path, int64_t start, int64_t end,
                 size_t buf_size = (1u << 20))
      : f_(std::fopen(path.c_str(), "rb")), buf_(buf_size), pos_(0), len_(0),
        remaining_(end < 0 ? std::numeric_limits<int64_t>::max() : end - start) {
    if (!f_) throw std::runtime_error("cannot open: " + path);
#if defined(__APPLE__)
    fcntl(fileno(f_), F_RDAHEAD, 1);
#elif defined(POSIX_FADV_SEQUENTIAL)
    posix_fadvise(fileno(f_), (off_t)start, end < 0 ? 0 : (off_t)(end - start),
                  POSIX_FADV_SEQUENTIAL);
#endif
    if (start > 0 && fseeko(f_, (off_t)start, SEEK_SET) != 0) {
      std::fclose(f_);
      f_ = nullptr;
      throw std::runtime_error("cannot seek in: " + path);
    }
  }

  ~FastLineReader() {
    if (f_) std::fclose(f_);
  }

  FastLineReader(const FastLineReader&) = delete;
  FastLineReader& operator=(const FastLineReader&) = delete;

  // read next line into out (without trailing '\n' / '\r')
  bool next_line(std::string& out) {
    out.clear();
    while (true) {
      if (pos_ >= len_) {
        size_t want = buf_.size();
        if ((int64_t)want > remaining_) want = (size_t)(remaining_ > 0 ? remaining_ : 0);
        len_ = want ? std::fread(buf_.data(), 1, want, f_) : 0;
        remaining_ -= (int64_t)len_;
        pos_ = 0;
        if (len_ == 0) return !out.empty();
      }
      const size_t start = pos_;
      const void* nl = std::memchr(buf_.data() + pos_, '\n', len_ - pos_);
      pos_ = nl ? (size_t)((const char*)nl - buf_.data()) : len_;
      out.append(buf_.data() + start, pos_ - start);
      if (nl) {
        pos_++;
        if (!out.empty() && out.back() == '\r') out.pop_back();
        return true;
      }
    }
  }

private:
  std::FILE* f_;
  std::vector<char> buf_;
  size_t pos_, len_;
  int64_t remaining_;
};

// merge key within one chromosome slice: byte-wise (kmerA_id, kmerB_id).
struct Key {
  std::string kA, kB;

  bool operator<(const Key& o) const {
    if (int c = kA.compare(o.kA); c != 0) return c < 0;
    return kB.compare(o.kB) < 0;
  }
  bool operator==(const Key& o) const { return kA == o.kA && kB == o.kB; }
};

struct KeyView {
  std::string_view kA, kB;

  bool operator<(const KeyView& o) const {
    if (int c = kA.compare(o.kA); c != 0) return c < 0;
    return kB.compare(o.kB) < 0;
  }
};

inline bool key_eq(const Key& k, const KeyView& v) {
  return std::string_view(k.kA) == v.kA && std::string_view(k.kB) == v.kB;
}

void split_tabs(std::string_view line, std::vector<std::string_view>& out) {
  out.clear();
  const char* base = line.data();
  const char* p = base;
  const char* end = base + line.size();
  while (true) {
    const char* t = (const char*)std::memchr(p, '\t', (size_t)(end - p));
    if (!t) break;
    out.push_back(std::string_view(p, (size_t)(t - p)));
    p = t + 1;
  }
  out.push_back(std::string_view(p, (size_t)(end - p)));
}

constexpr const char* DIGEST_HEADER =
    "chr\tkmerA_id\tkmerB_id\tn_reads\tallele_lengths\tmotif_kmers";
constexpr const char* DIGEST_SUFFIX = ".pon-digest.tsv";
constexpr const char* DETECT_SUFFIX = ".alleles.tsv";

enum class Format { ALLELE, DIGEST };

// column positions (0-based) of one input file, resolved from its header
struct Layout {
  Format fmt = Format::ALLELE;
  int chr = -1, kmer_a = -1, kmer_b = -1;
  int distance = -1, motif_count = -1;             // raw allele output (lists)
  int n_reads = -1, distances = -1, motif = -1;    // digest
  int n_cols = 0;  // a row with fewer fields than this is malformed and skipped
};

// validate a header line and resolve the layout.
Layout classify_header(std::string_view line, const std::string& path) {
  Layout L;
  int alen = -1;   // the allele_lengths column, whichever layout this is
  std::vector<std::string_view> f;
  split_tabs(line, f);
  for (int i = 0; i < (int)f.size(); i++) {
    std::string_view name = f[i];
    if (!name.empty() && name.back() == '\r') name.remove_suffix(1);
    if      (name == "chr")              L.chr = i;
    else if (name == "kmerA_id")         L.kmer_a = i;
    else if (name == "kmerB_id")         L.kmer_b = i;
    else if (name == "allele_lengths")   alen = i;
    else if (name == "n_motif_kmers")    L.motif_count = i;
    else if (name == "n_reads")          L.n_reads = i;
    else if (name == "motif_kmers")      L.motif = i;
  }
  const bool digest = L.motif >= 0;
  L.fmt = digest ? Format::DIGEST : Format::ALLELE;
  if (digest) L.distances = alen; else L.distance = alen;
  std::vector<std::pair<int, const char*>> need = {
      {L.chr, "chr"}, {L.kmer_a, "kmerA_id"}, {L.kmer_b, "kmerB_id"}};
  if (digest) {
    need.push_back({L.n_reads, "n_reads"});
    need.push_back({L.distances, "allele_lengths"});
    need.push_back({L.motif, "motif_kmers"});
  } else {
    need.push_back({L.distance, "allele_lengths"});
    need.push_back({L.motif_count, "n_motif_kmers"});
  }
  for (const auto& [idx, name] : need)
    if (idx < 0)
      throw std::runtime_error(path + ": header is missing the '" + name + "' column. Is this a "
                               "current trumut allele output or a build-pon digest?");
  L.n_cols = 1 + std::max({L.chr, L.kmer_a, L.kmer_b, L.distance, L.motif_count,
                           L.n_reads, L.distances, L.motif});
  return L;
}

bool parse_i32(std::string_view d, int32_t& out) {
  if (d.empty()) return false;
  size_t i = (d[0] == '-') ? 1 : 0;
  if (i >= d.size()) return false;
  int64_t val = 0;
  for (; i < d.size(); i++) {
    if (d[i] < '0' || d[i] > '9') return false;
    val = val * 10 + (d[i] - '0');
    if (val > std::numeric_limits<int32_t>::max()) return false;
  }
  out = (int32_t)(d[0] == '-' ? -val : val);
  return true;
}

bool parse_u64(std::string_view d, uint64_t& out) {
  if (d.empty()) return false;
  uint64_t val = 0;
  for (char c : d) {
    if (c < '0' || c > '9') return false;
    if (val > (std::numeric_limits<uint64_t>::max() - 9) / 10) return false;
    val = val * 10 + (uint64_t)(c - '0');
  }
  out = val;
  return true;
}

// split a comma-separated field into its items.
void split_commas(std::string_view s, std::vector<std::string_view>& out) {
  out.clear();
  const char* p = s.data();
  const char* end = p + s.size();
  while (true) {
    const char* c = (const char*)std::memchr(p, ',', (size_t)(end - p));
    if (!c) break;
    out.push_back(std::string_view(p, (size_t)(c - p)));
    p = c + 1;
  }
  out.push_back(std::string_view(p, (size_t)(end - p)));
}

struct Rep {
  int32_t dist;
  double motif = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> vals;
};

struct ChrBlock {
  ChrKey chr;
  int64_t start, end;
};

class FileIndexer {
public:
  explicit FileIndexer(const std::string& path) : path_(path) {
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) throw std::runtime_error("cannot open: " + path);
    struct stat st;
    if (fstat(fd_, &st) != 0) {
      ::close(fd_);
      throw std::runtime_error("cannot stat: " + path);
    }
    size_ = (int64_t)st.st_size;
  }

  ~FileIndexer() {
    if (fd_ >= 0) ::close(fd_);
  }

  FileIndexer(const FileIndexer&) = delete;
  FileIndexer& operator=(const FileIndexer&) = delete;

  int64_t read_and_check_header() {
    if (size_ == 0) throw std::runtime_error(path_ + ": empty file");
    lay_ = classify_header(read_line_at(0), path_);
    return next_line_start(1);  // just past the header's '\n'
  }
  const Layout& layout() const { return lay_; }

  std::vector<ChrBlock> find_blocks(int64_t data_start) {
    std::vector<ChrBlock> blocks;
    int64_t pos = data_start;
    while (pos < size_) {
      ChrKey k = chr_at(pos);
      if (!blocks.empty() && !(blocks.back().chr < k))
        throw std::runtime_error(path_ + ": chromosome blocks are not in karyotype "
                                 "order (" + blocks.back().chr.name + " then " + k.name +
                                 "). This tool requires the unmodified output of "
                                 "trumut allele on a BAM/CRAM input.");
      int64_t lo = pos + 1, hi = size_;
      while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        int64_t s = next_line_start(mid);
        if (s >= size_ || k < chr_at(s)) hi = mid;
        else lo = mid + 1;
      }
      int64_t end = next_line_start(lo);
      blocks.push_back({std::move(k), pos, end});
      pos = end;
    }
    return blocks;
  }

  bool kmerA_at_or_after(int64_t probe, int64_t end, std::string& out) {
    for (int64_t s = next_line_start(probe); s < end; s = next_line_start(s + 1))
      if (row_kmerA(s, out)) return true;
    return false;
  }

  int64_t find_split(int64_t start, int64_t end, const std::string& threshold) {
    std::string kA;
    auto at_or_past = [&](int64_t probe) {
      if (!kmerA_at_or_after(probe, end, kA)) return true;  // only tail junk left
      return !(kA < threshold);
    };
    int64_t lo = start, hi = end;
    while (lo < hi) {
      int64_t mid = lo + (hi - lo) / 2;
      if (at_or_past(mid)) hi = mid;
      else lo = mid + 1;
    }
    int64_t s = next_line_start(lo);
    return s < end ? s : end;
  }

private:

  int64_t next_line_start(int64_t off) {
    if (off <= 0) return 0;
    int64_t p = off - 1;  // a '\n' at off-1 makes off itself a line start
    char buf[1 << 16];
    while (p < size_) {
      ssize_t n = pread(fd_, buf, sizeof buf, (off_t)p);
      if (n < 0) throw std::runtime_error(path_ + ": read error");
      if (n == 0) break;
      for (ssize_t i = 0; i < n; i++)
        if (buf[i] == '\n') return p + i + 1;
      p += n;
    }
    return size_;
  }

  ChrKey chr_at(int64_t s) {
    char buf[4096];
    ssize_t n = pread(fd_, buf, sizeof buf, (off_t)s);
    if (n < 0) throw std::runtime_error(path_ + ": read error");
    ChrKey k;
    ssize_t i = 0;
    int col = 0;
    while (i < n && buf[i] != '\n' && col < lay_.chr) {
      if (buf[i] == '\t') col++;
      i++;
    }
    if (col == lay_.chr && i < n) {
      ssize_t j = i;
      while (j < n && buf[j] != '\t' && buf[j] != '\n') j++;
      if (j < n) {
        k.name.assign(buf + i, (size_t)(j - i));
        k.rank = chr_rank(k.name);
        return k;
      }
    }
    std::string line = read_line_at(s);
    split_tabs(line, probe_fields_);
    if ((int)probe_fields_.size() > lay_.chr) k.name.assign(probe_fields_[lay_.chr]);
    k.rank = chr_rank(k.name);
    return k;
  }

  bool row_kmerA(int64_t s, std::string& out) {
    char buf[4096];
    ssize_t n = pread(fd_, buf, sizeof buf, (off_t)s);
    if (n < 0) throw std::runtime_error(path_ + ": read error");
    if (n == 0) return false;
    const char* nl = (const char*)std::memchr(buf, '\n', (size_t)n);
    std::string spill;
    std::string_view line;
    if (nl) {
      line = std::string_view(buf, (size_t)(nl - buf));
    } else {
      spill = read_line_at(s);
      line = spill;
    }
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    split_tabs(line, probe_fields_);
    if ((int)probe_fields_.size() < lay_.n_cols) return false;
    out.assign(probe_fields_[lay_.kmer_a]);
    return true;
  }
  std::string read_line_at(int64_t s) {
    std::string out;
    char buf[1 << 16];
    int64_t p = s;
    while (p < size_) {
      ssize_t n = pread(fd_, buf, sizeof buf, (off_t)p);
      if (n < 0) throw std::runtime_error(path_ + ": read error");
      if (n == 0) break;
      for (ssize_t i = 0; i < n; i++) {
        if (buf[i] == '\n') {
          out.append(buf, (size_t)i);
          if (!out.empty() && out.back() == '\r') out.pop_back();
          return out;
        }
      }
      out.append(buf, (size_t)n);
      p += n;
      if (out.size() > (64u << 20)) throw std::runtime_error(path_ + ": header line too long");
    }
    if (!out.empty() && out.back() == '\r') out.pop_back();
    return out;
  }

  std::string path_;
  int fd_ = -1;
  int64_t size_ = 0;
  Layout lay_;
  std::vector<std::string_view> probe_fields_;
};

struct Source {
  FastLineReader rdr;
  const std::string& path;
  const std::string& expected_chr;
  const Layout lay;
  const Format fmt;
  Key key;                       // current group's key
  size_t n_reads = 0;
  std::vector<Rep> reps;
  KeyView pending_key;
  std::vector<int32_t> pending_lens;
  std::vector<double> pending_motifs;
  size_t pending_n_reads = 0;
  std::vector<Rep> pending_reps;
  bool has_pending = false;
  uint64_t bad_rows = 0;

  Source(const std::string& p, int64_t start, int64_t end, const std::string& chr,
         const Layout& layout, const std::string& bound_lo,
         const std::string& bound_hi)
      : rdr(p, start, end), path(p), expected_chr(chr), lay(layout), fmt(layout.fmt),
        bound_lo_(bound_lo), bound_hi_(bound_hi) {}

  bool next_row() {
    std::string& line = *scratch_;
    std::vector<std::string_view>& f = fields_scratch_;
    const bool digest = fmt == Format::DIGEST;
    const int n_cols = lay.n_cols;
    const int col_chr = lay.chr;
    const int col_kA = lay.kmer_a;
    const int col_kB = lay.kmer_b;
    while (rdr.next_line(line)) {
      split_tabs(line, f);
      if ((int)f.size() < n_cols) { bad_rows++; continue; }

      if (f[col_chr] != expected_chr)
        throw std::runtime_error(path + ": found a " + std::string(f[col_chr]) +
                                 " row inside the " + expected_chr + " block: rows are "
                                 "not in allele's alignment-input sort order. This tool "
                                 "requires the unmodified output of trumut allele on a "
                                 "BAM/CRAM input, or a build-pon digest of one.");

      if (digest) {
        if (!parse_digest_row(f)) { bad_rows++; continue; }
      } else {
        if (!parse_raw_lists(f)) { bad_rows++; continue; }
      }

      KeyView k{f[col_kA], f[col_kB]};
      if (!has_prev && !bound_lo_.empty() && k.kA < std::string_view(bound_lo_))
        throw std::runtime_error(path + ": " + expected_chr + " row " + std::string(k.kA) +
                                 " sorts before the start of its own kmerA_id slice (" +
                                 bound_lo_ + "): rows are not in allele's sort order. This "
                                 "tool requires the unmodified output of trumut allele.");
      if (has_prev && k < pending_key)
        throw std::runtime_error(path + ": " + expected_chr + " rows are not in "
                                 "allele's sort order (kmerA_id, kmerB_id). This tool "
                                 "requires the unmodified output of trumut allele.");
      has_prev = true;

      pending_key = k;
      has_pending = true;
      std::swap(scratch_, accepted_);
      return true;
    }
    if (has_prev && !bound_hi_.empty() &&
        !(pending_key.kA < std::string_view(bound_hi_)))
      throw std::runtime_error(path + ": " + expected_chr + " row " +
                               std::string(pending_key.kA) + " sorts at or past the end of "
                               "its own kmerA_id slice (" + bound_hi_ + "): rows are not in "
                               "allele's sort order. This tool requires the unmodified "
                               "output of trumut allele.");
    return false;
  }

  bool advance() {
    if (!has_pending && !next_row()) return false;
    key.kA.assign(pending_key.kA);
    key.kB.assign(pending_key.kB);
    reps.clear();

    if (fmt == Format::DIGEST) {
      n_reads = pending_n_reads;
      reps.swap(pending_reps);
      has_pending = false;
      if (next_row() && key_eq(key, pending_key))
        throw std::runtime_error(path + ": " + expected_chr + " pair (" + key.kA + ", " +
                                 key.kB + ") appears twice: a digest holds one row per "
                                 "pair. Is this file a valid build-pon digest?");
      return true;
    }

    n_reads = 0;
    auto take = [&]() {
      for (size_t i = 0; i < pending_lens.size(); i++) {
        n_reads++;
        Rep* r = nullptr;
        for (auto& x : reps)
          if (x.dist == pending_lens[i]) { r = &x; break; }
        if (!r) { reps.push_back(Rep{pending_lens[i]}); r = &reps.back(); }
        r->vals.push_back(pending_motifs[i]);
      }
      has_pending = false;
    };
    take();
    while (next_row()) {
      if (!key_eq(key, pending_key)) break;
      take();
    }
    for (auto& r : reps) { r.motif = median_of(r.vals); r.vals.clear(); }
    return true;
  }

private:
  bool parse_raw_lists(const std::vector<std::string_view>& f) {
    split_commas(f[lay.distance], items_a_);
    split_commas(f[lay.motif_count], items_b_);
    if (items_a_.empty()) return false;
    pending_lens.clear();
    pending_motifs.clear();
    for (size_t i = 0; i < items_a_.size(); i++) {
      int32_t d;
      if (!parse_i32(items_a_[i], d)) return false;
      pending_lens.push_back(d);
      pending_motifs.push_back(i < items_b_.size() ? parse_motif_count(items_b_[i])
                                                    : std::numeric_limits<double>::quiet_NaN());
    }
    return true;
  }

  bool parse_digest_row(const std::vector<std::string_view>& f) {
    uint64_t n;
    if (!parse_u64(f[lay.n_reads], n)) return false;
    split_commas(f[lay.distances], items_a_);
    const bool na = f[lay.motif] == "NA";
    if (!na) {
      split_commas(f[lay.motif], items_b_);
      if (items_b_.size() != items_a_.size()) return false;
    }
    pending_reps.clear();
    pending_reps.reserve(items_a_.size());
    for (size_t i = 0; i < items_a_.size(); i++) {
      Rep r{0};
      if (!parse_i32(items_a_[i], r.dist)) return false;
      if (!na) {
        r.motif = parse_motif_count(items_b_[i]);
        if (std::isnan(r.motif) && items_b_[i] != "NA") return false;
      }
      for (const auto& prev : pending_reps)
        if (prev.dist == r.dist) return false;
      pending_reps.push_back(std::move(r));
    }
    if (pending_reps.empty()) return false;
    pending_n_reads = (size_t)std::min<uint64_t>(n, std::numeric_limits<size_t>::max());
    return true;
  }

  bool has_prev = false;
  const std::string& bound_lo_;
  const std::string& bound_hi_;
  std::string line_a_, line_b_;
  std::string* scratch_ = &line_a_;
  std::string* accepted_ = &line_b_;
  std::vector<std::string_view> fields_scratch_;
  std::vector<std::string_view> items_a_, items_b_;  // digest list-field scratch
};

struct HeapCmp {
  const std::vector<std::unique_ptr<Source>>* srcs;
  bool operator()(int a, int b) const {
    return (*srcs)[b]->key < (*srcs)[a]->key;
  }
};

struct ChrTask {
  ChrKey chr;
  std::string bound_lo, bound_hi;
  std::vector<int> file_idx;                       // files with rows in this slice
  std::vector<std::pair<int64_t, int64_t>> range;  // parallel to file_idx
  // results
  std::string out;                                 // finished output rows
  uint64_t rows_out = 0, rows_dropped = 0;
  std::vector<std::pair<int, uint64_t>> bad_rows;  // (file index, count)
  std::string error;                               // non-empty = task failed
};

void merge_chromosome(ChrTask& task, const std::vector<std::string>& paths,
                      const std::vector<Layout>& layouts, int n_files_total,
                      double pair_percent_support, double distance_percent_support,
                      long min_reads) {
  const int M = (int)task.file_idx.size();
  std::vector<std::unique_ptr<Source>> srcs;
  srcs.reserve(M);
  for (int i = 0; i < M; i++)
    srcs.push_back(std::make_unique<Source>(paths[task.file_idx[i]],
                                            task.range[i].first, task.range[i].second,
                                            task.chr.name, layouts[task.file_idx[i]],
                                            task.bound_lo, task.bound_hi));

  task.out.reserve(64u << 10);

  HeapCmp cmp{&srcs};
  std::priority_queue<int, std::vector<int>, HeapCmp> pq(cmp);
  for (int i = 0; i < M; i++)
    if (srcs[i]->advance()) pq.push(i);

  const auto meets = [](size_t count, size_t total, double frac) {
    return (double)count + 1e-9 >= frac * (double)total;
  };

  std::vector<int32_t> pair_dists;
  pair_dists.reserve(256);
  size_t covering = 0;  // covering files of the current pair
  std::vector<Rep> pair_reps;
  pair_reps.reserve(256);
  std::vector<double> med_counts;  // per-file motif medians being medianed

  Key current;
  bool have_key = false;

  auto append_i32 = [](std::string& s, int32_t v) {
    char b[16];
    auto res = std::to_chars(b, b + sizeof b, v);
    s.append(b, (size_t)(res.ptr - b));
  };

  // median across the covering files of their per-file median n_motif_kmers at distance D
  auto append_median_at = [&](std::string& o, int32_t D) {
    med_counts.clear();
    for (const auto& r : pair_reps)
      if (r.dist == D) med_counts.push_back(r.motif);
    o += median_string(med_counts, "NA");
  };

  auto flush_pair = [&]() {
    const size_t n_cov = covering;
    covering = 0;
    if (!meets(n_cov, (size_t)n_files_total, pair_percent_support)) {
      task.rows_dropped++;
      pair_dists.clear();
      pair_reps.clear();
      return;
    }

    std::sort(pair_dists.begin(), pair_dists.end());
    int32_t smin = 0, smax = 0;
    bool have_supported = false;
    for (size_t i = 0, j; i < pair_dists.size(); i = j) {
      for (j = i + 1; j < pair_dists.size() && pair_dists[j] == pair_dists[i]; j++) {}
      if (meets(j - i, n_cov, distance_percent_support)) {  // run length = files carrying this distance
        if (!have_supported) smin = pair_dists[i];
        smax = pair_dists[i];
        have_supported = true;
      }
    }
    pair_dists.clear();

    std::string& o = task.out;
    o += task.chr.name;
    o += '\t';
    o += current.kA;
    o += '\t';
    o += current.kB;
    if (have_supported) {
      o += '\t';
      append_i32(o, smin);
      o += '\t';
      append_i32(o, smax);
      o += '\t';
      append_median_at(o, smin);
      o += '\t';
      append_median_at(o, smax);
    } else {
      o += "\tNA\tNA\tNA\tNA";
    }
    o += '\n';
    task.rows_out++;
    pair_reps.clear();
  };

  while (!pq.empty()) {
    int si = pq.top();
    pq.pop();
    Source& s = *srcs[si];

    if (!have_key) {
      current = s.key;
      have_key = true;
    } else if (!(s.key == current)) {
      flush_pair();
      current = s.key;
    }

    if (s.n_reads >= (size_t)min_reads) {
      for (auto& r : s.reps) {
        pair_dists.push_back(r.dist);
        pair_reps.push_back(std::move(r));
      }
      covering++;
    }

    if (s.advance()) pq.push(si);
  }
  if (have_key) flush_pair();

  for (int i = 0; i < M; i++)
    if (srcs[i]->bad_rows) task.bad_rows.emplace_back(task.file_idx[i], srcs[i]->bad_rows);
}

struct DigestStats {
  uint64_t rows_in = 0, pairs_out = 0, bad_rows = 0;
};

// reduce one raw allele output to its digest.
DigestStats digest_file(const std::string& in, const std::string& out) {
  DigestStats st;
  FastLineReader rdr(in, 0, -1, 8u << 20);
  std::string line;
  if (!rdr.next_line(line)) throw std::runtime_error(in + ": empty file");
  const Layout lay = classify_header(line, in);
  if (lay.fmt != Format::ALLELE)
    throw std::runtime_error(in + ": is already a digest");

  const std::string tmp = out + ".tmp";
  std::FILE* fo = std::fopen(tmp.c_str(), "wb");
  if (!fo) throw std::runtime_error("cannot open output: " + tmp);
  std::vector<char> obuf(4u << 20);
  std::setvbuf(fo, obuf.data(), _IOFBF, obuf.size());
  auto fail_write = [&]() {
    std::fclose(fo);
    std::remove(tmp.c_str());
    throw std::runtime_error("write error: " + tmp);
  };
  if (std::fputs(DIGEST_HEADER, fo) < 0 || std::fputc('\n', fo) == EOF) fail_write();

  std::vector<std::string_view> f;
  std::string row;

  // current chromosome block
  ChrKey cur_chr;
  bool have_chr = false;

  // current group
  Key key;
  bool in_group = false;
  size_t n_reads = 0;            // reads in the group (list items)
  std::vector<Rep> reps;         // every read's n_motif_kmers at each distinct distance
  std::vector<std::string_view> lens_items, motif_items;
  std::vector<int32_t> row_lens;

  auto append_u64 = [](std::string& s, uint64_t v) {
    char b[24];
    auto res = std::to_chars(b, b + sizeof b, v);
    s.append(b, (size_t)(res.ptr - b));
  };
  auto append_i32 = [](std::string& s, int32_t v) {
    char b[16];
    auto res = std::to_chars(b, b + sizeof b, v);
    s.append(b, (size_t)(res.ptr - b));
  };

  auto flush_group = [&]() {
    if (!in_group) return;

    row.clear();
    row += cur_chr.name;
    row += '\t';
    row += key.kA;
    row += '\t';
    row += key.kB;
    row += '\t';
    append_u64(row, n_reads);
    row += '\t';
    for (size_t i = 0; i < reps.size(); i++) {
      if (i) row += ',';
      append_i32(row, reps[i].dist);
    }
    row += '\t';
    // per-distance median of the reads' n_motif_kmers
    bool all_na = true;
    for (auto& r : reps) { r.motif = median_of(r.vals); if (!std::isnan(r.motif)) all_na = false; }
    if (all_na) {
      row += "NA";
    } else {
      for (size_t i = 0; i < reps.size(); i++) {
        if (i) row += ',';
        row += format_median(reps[i].motif, "NA");
      }
    }
    row += '\n';
    if (std::fwrite(row.data(), 1, row.size(), fo) != row.size()) fail_write();
    st.pairs_out++;
    in_group = false;
    n_reads = 0;
    reps.clear();
  };

  auto order_error = [&](const std::string& what) {
    std::fclose(fo);
    std::remove(tmp.c_str());
    return std::runtime_error(in + ": " + what + ". This tool requires the unmodified "
                              "output of trumut allele on a BAM/CRAM input.");
  };

  while (rdr.next_line(line)) {
    st.rows_in++;
    split_tabs(line, f);
    if ((int)f.size() < lay.n_cols) { st.bad_rows++; continue; }
    split_commas(f[lay.distance], lens_items);
    split_commas(f[lay.motif_count], motif_items);
    if (lens_items.empty()) { st.bad_rows++; continue; }
    row_lens.clear();
    bool bad = false;
    for (auto it : lens_items) { int32_t d; if (!parse_i32(it, d)) { bad = true; break; } row_lens.push_back(d); }
    if (bad) { st.bad_rows++; continue; }

    if (!have_chr || f[lay.chr] != cur_chr.name) {
      flush_group();
      ChrKey k;
      k.name.assign(f[lay.chr]);
      k.rank = chr_rank(k.name);
      if (have_chr && !(cur_chr < k))
        throw order_error("chromosome blocks are not in karyotype order (" + cur_chr.name +
                          " then " + k.name + ")");
      cur_chr = std::move(k);
      have_chr = true;
    }

    KeyView kv{f[lay.kmer_a], f[lay.kmer_b]};
    if (in_group && kv < KeyView{key.kA, key.kB})
      throw order_error(cur_chr.name + " rows are not in allele's sort order (kmerA_id, kmerB_id)");

    if (!in_group || !key_eq(key, kv)) {
      flush_group();
      key.kA.assign(kv.kA);
      key.kB.assign(kv.kB);
      in_group = true;
    }
    for (size_t i = 0; i < row_lens.size(); i++) {
      n_reads++;
      Rep* r = nullptr;
      for (auto& x : reps)
        if (x.dist == row_lens[i]) { r = &x; break; }
      if (!r) { reps.push_back(Rep{row_lens[i]}); r = &reps.back(); }
      r->vals.push_back(i < motif_items.size() ? parse_motif_count(motif_items[i])
                                               : std::numeric_limits<double>::quiet_NaN());
    }
  }
  flush_group();

  if (std::fclose(fo) != 0) {
    std::remove(tmp.c_str());
    throw std::runtime_error("write error: " + tmp);
  }
  if (std::rename(tmp.c_str(), out.c_str()) != 0) {
    std::remove(tmp.c_str());
    throw std::runtime_error("cannot rename " + tmp + " to " + out);
  }
  return st;
}

// sample name: filename minus the given suffix.
std::string stem_of(const fs::path& p, const std::string& suffix) {
  std::string name = p.filename().string();
  return name.substr(0, name.size() - suffix.size());
}

uint64_t raise_nofile_limit(uint64_t want) {
  struct rlimit rl;
  if (getrlimit(RLIMIT_NOFILE, &rl) != 0) return 1024;
  uint64_t hard = (rl.rlim_max == RLIM_INFINITY) ? want : (uint64_t)rl.rlim_max;
#ifdef __APPLE__
  if (hard > 10240) hard = 10240;
#endif
  uint64_t target = want < hard ? want : hard;
  if ((uint64_t)rl.rlim_cur < target) {
    rl.rlim_cur = (rlim_t)target;
    if (setrlimit(RLIMIT_NOFILE, &rl) != 0) {
      getrlimit(RLIMIT_NOFILE, &rl);
      return (uint64_t)rl.rlim_cur;
    }
  }
  getrlimit(RLIMIT_NOFILE, &rl);
  return (uint64_t)rl.rlim_cur;
}

int digest_main(const std::vector<const char*>& pos, const std::string& suffix, long jobs) {
  if (pos.empty() || pos.size() > 2) {
    throw std::runtime_error("--digest takes <in.alleles.tsv> <out.pon-digest.tsv>, "
                             "or <dir> [<out_dir>]");
  }
  const fs::path in = pos[0];
  if (!fs::exists(in)) throw std::runtime_error(std::string(pos[0]) + ": no such file or directory");

  struct Job { std::string in, out; };
  std::vector<Job> work;
  if (fs::is_regular_file(in)) {
    if (pos.size() != 2)
      throw std::runtime_error("--digest <file> needs an output path");
    work.push_back({in.string(), pos[1]});
  } else if (fs::is_directory(in)) {
    const fs::path out_dir = pos.size() == 2 ? fs::path(pos[1]) : in;
    if (!fs::exists(out_dir) || !fs::is_directory(out_dir))
      throw std::runtime_error(out_dir.string() + " is not a directory");
    std::vector<fs::path> files;
    for (auto& e : fs::directory_iterator(in))
      if (e.is_regular_file() && ends_with(e.path().filename().string(), suffix))
        files.push_back(e.path());
    std::sort(files.begin(), files.end());
    if (files.empty()) throw std::runtime_error("no *" + suffix + " files found in " + in.string());
    int skipped = 0;
    for (auto& f : files) {
      const fs::path out = out_dir / (stem_of(f, suffix) + DIGEST_SUFFIX);
      std::error_code ec;
      if (fs::exists(out, ec) && fs::last_write_time(out, ec) >= fs::last_write_time(f, ec) && !ec) {
        skipped++;
        continue;
      }
      work.push_back({f.string(), out.string()});
    }
    std::cerr << "[INFO] found " << files.size() << " *" << suffix << " file(s); " << skipped
              << " already digested, " << work.size() << " to do\n";
    if (work.empty()) return 0;
  } else {
    throw std::runtime_error(std::string(pos[0]) + ": not a regular file or directory");
  }

  long J = jobs > 0 ? jobs : 3;
  if (J > (long)work.size()) J = (long)work.size();

  std::cerr << "[INFO] digesting " << J << " file(s) at a time\n";
  const auto t0 = std::chrono::steady_clock::now();
  auto elapsed = [&]() {
    long s = (long)std::chrono::duration_cast<std::chrono::seconds>(
                 std::chrono::steady_clock::now() - t0).count();
    char b[32];
    std::snprintf(b, sizeof b, "%ld:%02ld:%02ld", s / 3600, (s / 60) % 60, s % 60);
    return std::string(b);
  };

  std::mutex log_mu;
  std::atomic<int> next{0};
  std::atomic<int> failed{0};
  std::atomic<int> done{0};
  std::vector<std::thread> pool;
  for (long t = 0; t < J; t++) {
    pool.emplace_back([&]() {
      for (;;) {
        int i = next.fetch_add(1);
        if (i >= (int)work.size()) return;
        {
          std::lock_guard<std::mutex> lk(log_mu);
          std::cerr << "[INFO] " << elapsed() << " start " << work[i].in << "\n";
        }
        try {
          DigestStats st = digest_file(work[i].in, work[i].out);
          const int k = ++done;
          std::lock_guard<std::mutex> lk(log_mu);
          std::cerr << "[INFO] " << elapsed() << " done " << k << "/" << work.size() << " "
                    << work[i].in << ": " << st.rows_in << " rows -> " << st.pairs_out
                    << " pairs";
          if (st.bad_rows) std::cerr << " (skipped " << st.bad_rows << " malformed rows)";
          std::cerr << "\n";
        } catch (const std::exception& e) {
          failed++;
          ++done;
          std::lock_guard<std::mutex> lk(log_mu);
          std::cerr << "[ERROR] " << work[i].in << ": " << e.what() << "\n";
        }
      }
    });
  }
  for (auto& th : pool) th.join();
  std::cerr << "[INFO] finished " << work.size() << " file(s) in " << elapsed();
  if (failed) {
    std::cerr << ", " << failed.load() << " failed\n";
    return 1;
  }
  std::cerr << "\n";
  return 0;
}

void usage(std::FILE* fp, const char* prog) {
  (void)prog;
  std::fprintf(fp,
        "Usage: trumut build-pon alleles_dir/ pon.tsv --input-type reads|assembly [options]\n"
        "       trumut build-pon --digest sample.alleles.tsv sample.pon-digest.tsv\n"
        "       trumut build-pon --digest alleles_dir/ [out_dir/]\n"
        "\n"
        "Positional:\n"
        "  alleles_dir/     directory holding one trumut allele output per normal sample, named *.alleles.tsv\n"
        "  pon.tsv          output path of the panel of normals\n"
        "\n"
        "Required (merge):\n"
        "  --input-type STR              reads for sequencing reads, assembly for assembled contigs\n"
        "  --ploidy N                    required with assembly: haplotype contigs per locus (1 haploid, 2 diploid)\n"
        "\n"
        "Options:\n"
        "  --pair-percent-support F      fraction of ALL samples that must cover a pair for it to be kept [0.98]\n"
        "  --distance-percent-support F  fraction of the COVERING samples that must carry a length for it to count as supported [0.0975]\n"
        "  --min-reads-per-sample N      a sample covers a pair only with at least N reads / contigs there [10; for assemblies, set by --ploidy]\n"
        "  --suffix S                    input filename suffix [.alleles.tsv]\n"
        "  --jobs J                      merge tasks at once, or files digested at once in --digest directory mode [3]\n"
        "  --split-factor N              cut every chromosome into exactly N kmerA_id sub-slices for parallelization [by chromosome size]\n"
        "  --digest                      reduce raw outputs to smaller digests that merge to the same panel; the merge\n"
        "                                uses a sample's digest when one exists next to its allele output\n"
        "  -h, --help                    show this message\n"
        "\n"
        "Output columns: chr kmerA_id kmerB_id range_min range_max range_min_n_motif_kmers range_max_n_motif_kmers\n"
        "\n"
        "Documentation: https://erkim11.github.io/trumut/\n");
}

}  // namespace

int build_pon_main(int argc, char** argv) {
  try {
    if (argc > 1 && (std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0)) {
      usage(stdout, argv[0]);
      return 0;
    }

    std::vector<const char*> pos;
    double pair_percent_support = 0.98;
    double distance_percent_support = 0.0975;
    std::string suffix = DETECT_SUFFIX;
    bool suffix_set = false;
    bool digest_mode = false;
    long min_reads = -1;   // --min-reads-per-sample; -1 = not given
    int input_type = 0;    // --input-type: 1 reads, -1 assembly; required for the merge
    int ploidy = 0;        // --ploidy: 0 = not given
    long jobs = 0;  // 0 = auto
    long split_factor = 0;  // 0 = auto (size-based); >0 forces sub-slices per chromosome
    for (int i = 1; i < argc; i++) {
      auto need_value = [&](const char* flag) -> const char* {
        if (i + 1 >= argc)
          throw std::runtime_error(std::string(flag) + " requires a value");
        return argv[++i];
      };
      if (std::strcmp(argv[i], "--pair-percent-support") == 0) {
        char* end = nullptr;
        pair_percent_support = std::strtod(need_value("--pair-percent-support"), &end);
        if (!end || *end || !(pair_percent_support > 0.0) || pair_percent_support > 1.0)
          throw std::runtime_error("--pair-percent-support must be a number in (0, 1]");
      } else if (std::strcmp(argv[i], "--distance-percent-support") == 0) {
        char* end = nullptr;
        distance_percent_support = std::strtod(need_value("--distance-percent-support"), &end);
        if (!end || *end || !(distance_percent_support > 0.0) || distance_percent_support > 1.0)
          throw std::runtime_error("--distance-percent-support must be a number in (0, 1]");
      } else if (std::strcmp(argv[i], "--suffix") == 0) {
        suffix = need_value("--suffix");
        suffix_set = true;
      } else if (std::strcmp(argv[i], "--digest") == 0) {
        digest_mode = true;
      } else if (std::strcmp(argv[i], "--min-reads-per-sample") == 0) {
        char* end = nullptr;
        min_reads = std::strtol(need_value("--min-reads-per-sample"), &end, 10);
        if (!end || *end || min_reads < 1)
          throw std::runtime_error("--min-reads-per-sample must be an integer >= 1");
      } else if (std::strcmp(argv[i], "--input-type") == 0) {
        const char* v = need_value("--input-type");
        input_type = parse_input_type(v);
        if (input_type == 0)
          throw std::runtime_error(std::string("--input-type must be 'reads' or 'assembly' (got '") + v + "')");
      } else if (std::strcmp(argv[i], "--ploidy") == 0) {
        char* end = nullptr;
        ploidy = (int)std::strtol(need_value("--ploidy"), &end, 10);
        if (!end || *end || ploidy < 1)
          throw std::runtime_error("--ploidy must be a positive integer (1 haploid, 2 diploid, ...)");
      } else if (std::strcmp(argv[i], "--jobs") == 0 || std::strcmp(argv[i], "-j") == 0) {
        char* end = nullptr;
        jobs = std::strtol(need_value("--jobs"), &end, 10);
        if (!end || *end || jobs < 1)
          throw std::runtime_error("--jobs must be an integer >= 1");
      } else if (std::strcmp(argv[i], "--split-factor") == 0) {
        char* end = nullptr;
        split_factor = std::strtol(need_value("--split-factor"), &end, 10);
        if (!end || *end || split_factor < 1)
          throw std::runtime_error("--split-factor must be an integer >= 1");
      } else if (argv[i][0] == '-' && argv[i][1]) {
        throw std::runtime_error(std::string("unknown option: ") + argv[i]);
      } else {
        pos.push_back(argv[i]);
      }
    }
    if (digest_mode) return digest_main(pos, suffix, jobs);

    if (pos.size() != 2) {
      usage(stderr, argv[0]);
      return 2;
    }
    if (input_type == 0) {
      std::cerr << "[ERROR] " << INPUT_TYPE_REQUIRED_MSG;
      return 2;
    }
    const bool assembly = input_type == -1;
    {
      const size_t floor = resolve_min_reads(assembly, ploidy, (int)min_reads, "[ERROR] ");
      if (floor == MIN_READS_INVALID) return 2;
      min_reads = (long)floor;
    }
    std::cerr << "[INFO] panel of " << (assembly ? "assembly" : "read") << " samples; a file covers a pair with at least "
              << min_reads << (assembly ? " haplotype contigs" : " reads") << "\n";
    const fs::path input_dir = pos[0];
    const std::string out_tsv = pos[1];

    if (!fs::exists(input_dir) || !fs::is_directory(input_dir))
      throw std::runtime_error("input_dir is not a directory");

    std::map<std::string, fs::path> by_stem;
    int n_digest = 0;
    for (auto& e : fs::directory_iterator(input_dir)) {
      if (!e.is_regular_file()) continue;
      const std::string name = e.path().filename().string();
      if (suffix_set) {
        if (ends_with(name, suffix)) by_stem.emplace(name, e.path());
        continue;
      }
      if (ends_with(name, DIGEST_SUFFIX)) {
        by_stem[stem_of(e.path(), DIGEST_SUFFIX)] = e.path();  // overrides a raw file
      } else if (ends_with(name, DETECT_SUFFIX)) {
        by_stem.emplace(stem_of(e.path(), DETECT_SUFFIX), e.path());  // keeps a digest
      }
    }
    std::vector<fs::path> files;
    for (auto& [stem, p] : by_stem) {
      files.push_back(p);
      if (ends_with(p.filename().string(), DIGEST_SUFFIX)) n_digest++;
    }
    std::sort(files.begin(), files.end());
    if (files.empty())
      throw std::runtime_error(suffix_set ? "no *" + suffix + " files found in input_dir"
                                          : std::string("no *") + DETECT_SUFFIX + " or *" +
                                                DIGEST_SUFFIX + " files found in input_dir");

    const int N = (int)files.size();
    std::vector<std::string> paths;
    paths.reserve(N);
    for (auto& f : files) paths.push_back(f.string());
    std::cerr << "[INFO] found " << N << " sample(s)";
    if (!suffix_set) std::cerr << " (" << n_digest << " digest(s), " << (N - n_digest) << " raw)";
    std::cerr << "\n";

    // index every file's chromosome blocks (binary search)
    std::map<ChrKey, std::vector<std::pair<int, std::pair<int64_t, int64_t>>>> by_chr;
    std::vector<std::unique_ptr<FileIndexer>> indexers;
    std::vector<Layout> layouts(N);
    indexers.reserve(N);
    for (int i = 0; i < N; i++) {
      indexers.push_back(std::make_unique<FileIndexer>(paths[i]));
      int64_t data_start = indexers[i]->read_and_check_header();
      layouts[i] = indexers[i]->layout();
      for (auto& b : indexers[i]->find_blocks(data_start))
        by_chr[b.chr].push_back({i, {b.start, b.end}});
    }
    const int n_chr = (int)by_chr.size();

    long J = jobs > 0 ? jobs : 3;

    int64_t grand_total = 0;
    for (auto& [chr, slices] : by_chr)
      for (auto& [fi, r] : slices) grand_total += r.second - r.first;
    const int64_t target = std::max<int64_t>(1, grand_total / std::max(1L, J * 4));

    std::vector<ChrTask> tasks;
    for (auto& [chr, slices] : by_chr) {
      int64_t chr_bytes = 0;
      size_t biggest = 0;
      for (size_t i = 0; i < slices.size(); i++) {
        const int64_t sz = slices[i].second.second - slices[i].second.first;
        chr_bytes += sz;
        if (sz > slices[biggest].second.second - slices[biggest].second.first) biggest = i;
      }
      int S = split_factor > 0 ? (int)split_factor
                               : (int)std::min<int64_t>(64, std::max<int64_t>(1, chr_bytes / target));

      // thresholds sampled from the largest file's block
      std::vector<std::string> thresholds;
      if (S > 1) {
        const int fi = slices[biggest].first;
        const int64_t bs = slices[biggest].second.first, be = slices[biggest].second.second;
        std::string kA;
        for (int k = 1; k < S; k++) {
          const int64_t probe = bs + (be - bs) * k / S;
          if (indexers[fi]->kmerA_at_or_after(probe, be, kA)) thresholds.push_back(kA);
        }
        std::sort(thresholds.begin(), thresholds.end());
        thresholds.erase(std::unique(thresholds.begin(), thresholds.end()), thresholds.end());
      }

      // cut every file's block at each threshold
      const size_t n_sub = thresholds.size() + 1;
      std::vector<ChrTask> sub(n_sub);
      for (size_t j = 0; j < n_sub; j++) {
        sub[j].chr = chr;
        if (j > 0) sub[j].bound_lo = thresholds[j - 1];
        if (j + 1 < n_sub) sub[j].bound_hi = thresholds[j];
      }
      for (auto& [fi, r] : slices) {
        int64_t prev = r.first;
        for (size_t j = 0; j < n_sub; j++) {
          const int64_t cut = (j + 1 < n_sub)
                                  ? indexers[fi]->find_split(prev, r.second, thresholds[j])
                                  : r.second;
          if (cut > prev) {
            sub[j].file_idx.push_back(fi);
            sub[j].range.push_back({prev, cut});
          }
          prev = cut;
        }
      }
      for (auto& t : sub)
        if (!t.file_idx.empty()) tasks.push_back(std::move(t));
    }
    indexers.clear();  // release the index FDs before the workers claim theirs
    const int n_tasks = (int)tasks.size();
    if (J > n_tasks) J = n_tasks;
    uint64_t fd_limit = raise_nofile_limit((uint64_t)J * (uint64_t)N + 64);
    long J_fd = (long)((fd_limit > 64 ? fd_limit - 64 : 1) / (uint64_t)N);
    if (J_fd < 1) J_fd = 1;
    if (J > J_fd) {
      std::cerr << "[INFO] open-file limit " << fd_limit << " caps parallelism at "
                << J_fd << " chromosome(s) at once (" << N << " files each)\n";
      J = J_fd;
    }
    std::cerr << "[INFO] merging " << n_chr << " chromosome(s) as " << n_tasks << " task(s), "
              << J << " at a time\n";

    // run the per-chromosome merges
    std::atomic<int> next_task{0};
    std::vector<std::thread> pool;
    pool.reserve((size_t)J);
    for (long t = 0; t < J; t++) {
      pool.emplace_back([&]() {
        for (;;) {
          int ti = next_task.fetch_add(1);
          if (ti >= n_tasks) return;
          try {
            merge_chromosome(tasks[ti], paths, layouts, N, pair_percent_support,
                             distance_percent_support, min_reads);
          } catch (const std::exception& e) {
            tasks[ti].error = e.what();
          }
        }
      });
    }
    for (auto& th : pool) th.join();

    for (const auto& t : tasks)
      if (!t.error.empty())
        throw std::runtime_error(t.chr.name + ": " + t.error);

    // write the results in karyotype order
    std::FILE* out = std::fopen(out_tsv.c_str(), "wb");
    if (!out) throw std::runtime_error("cannot open output: " + out_tsv);
    std::fputs("chr\tkmerA_id\tkmerB_id\trange_min\trange_max\t"
               "range_min_n_motif_kmers\trange_max_n_motif_kmers\n", out);
    uint64_t rows_out = 0, rows_dropped = 0;
    for (const auto& t : tasks) {
      if (!t.out.empty() && std::fwrite(t.out.data(), 1, t.out.size(), out) != t.out.size())
        throw std::runtime_error("write error: " + out_tsv);
      rows_out += t.rows_out;
      rows_dropped += t.rows_dropped;
    }
    if (std::fclose(out) != 0) throw std::runtime_error("write error: " + out_tsv);

    std::vector<uint64_t> bad_by_file(N, 0);
    for (const auto& t : tasks)
      for (const auto& [fi, n] : t.bad_rows) bad_by_file[fi] += n;
    uint64_t bad_total = 0;
    for (int i = 0; i < N; i++) {
      if (bad_by_file[i])
        std::cerr << "[WARNING] " << paths[i] << ": skipped " << bad_by_file[i]
                  << " malformed rows\n";
      bad_total += bad_by_file[i];
    }
    std::cerr << "[INFO] wrote " << rows_out << " pair rows; dropped " << rows_dropped
              << " covered by < " << pair_percent_support * 100.0
              << "% of files" << (bad_total ? " (see warnings above)" : "") << "\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "[ERROR] " << e.what() << "\n";
    return 1;
  }
}
