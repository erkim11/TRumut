/*
prep-catalog.cpp -- `trumut prep-catalog`: builds a uniqmer catalog for trumut allele from your own tandem repeat catalog.
*/

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <zlib.h>

#include "third_party/ketopt.h"

namespace fs = std::filesystem;

namespace {

constexpr int64_t KMER_LEN = 31;

constexpr const char* OUT_SUFFIX_HEADER =
    "kmerA_ref_coords\tkmerA_id\tkmerA_seq\tkmerA_strand\t"
    "kmerB_ref_coords\tkmerB_id\tkmerB_seq\tkmerB_strand";

struct KrefLayout {
  int chr = 0, pos = 1, id = 3, strand = 5;
  int n_cols = 6;
  bool has_header = false;
};

class LineSource {
public:
  virtual ~LineSource() = default;
  virtual bool next_line(std::string& out) = 0;
};

class GzLineReader : public LineSource {
public:
  explicit GzLineReader(const std::string& path, size_t buf_size = (1u << 20))
      : f_(gzopen(path.c_str(), "rb")), buf_(buf_size), pos_(0), len_(0) {
    if (!f_) throw std::runtime_error("cannot open: " + path);
  }

  ~GzLineReader() override {
    if (f_) gzclose(f_);
  }

  GzLineReader(const GzLineReader&) = delete;
  GzLineReader& operator=(const GzLineReader&) = delete;

  bool next_line(std::string& out) override {
    out.clear();
    while (true) {
      if (pos_ >= len_) {
        int n = gzread(f_, buf_.data(), (unsigned)buf_.size());
        if (n < 0) throw std::runtime_error("read error (bad gzip?)");
        len_ = (size_t)n;
        pos_ = 0;
        if (len_ == 0) return !out.empty();
      }
      size_t start = pos_;
      while (pos_ < len_ && buf_[pos_] != '\n') pos_++;
      out.append(buf_.data() + start, pos_ - start);
      if (pos_ < len_ && buf_[pos_] == '\n') {
        pos_++;
        if (!out.empty() && out.back() == '\r') out.pop_back();
        return true;
      }
    }
  }

private:
  gzFile f_;
  std::vector<char> buf_;
  size_t pos_, len_;
};

class RangeFileReader : public LineSource {
public:
  RangeFileReader(const std::string& path,
                  std::vector<std::pair<int64_t, int64_t>> ranges,
                  size_t buf_size = (1u << 20))
      : f_(std::fopen(path.c_str(), "rb")), ranges_(std::move(ranges)), buf_(buf_size) {
    if (!f_) throw std::runtime_error("cannot open: " + path);
  }

  ~RangeFileReader() override {
    if (f_) std::fclose(f_);
  }

  RangeFileReader(const RangeFileReader&) = delete;
  RangeFileReader& operator=(const RangeFileReader&) = delete;

  bool next_line(std::string& out) override {
    out.clear();
    while (true) {
      if (pos_ >= len_) {
        if (!refill()) return !out.empty();
      }
      size_t start = pos_;
      while (pos_ < len_ && buf_[pos_] != '\n') pos_++;
      out.append(buf_.data() + start, pos_ - start);
      if (pos_ < len_ && buf_[pos_] == '\n') {
        pos_++;
        if (!out.empty() && out.back() == '\r') out.pop_back();
        return true;
      }
    }
  }

private:
  std::FILE* f_;
  std::vector<std::pair<int64_t, int64_t>> ranges_;
  size_t next_range_ = 0;
  int64_t cur_ = 0, end_ = 0;  // byte cursor within the current range
  std::vector<char> buf_;
  size_t pos_ = 0, len_ = 0;

  bool refill() {
    while (cur_ >= end_) {
      if (next_range_ >= ranges_.size()) return false;
      cur_ = ranges_[next_range_].first;
      end_ = ranges_[next_range_].second;
      next_range_++;
      if (fseeko(f_, (off_t)cur_, SEEK_SET) != 0)
        throw std::runtime_error("seek failed in kmer_refs");
    }
    size_t want = (size_t)std::min<int64_t>((int64_t)buf_.size(), end_ - cur_);
    len_ = std::fread(buf_.data(), 1, want, f_);
    pos_ = 0;
    if (len_ == 0) throw std::runtime_error("kmer_refs shrank while indexed ranges were being read");
    cur_ += (int64_t)len_;
    return true;
  }
};

void split_tabs(std::string_view line, std::vector<std::string_view>& out) {
  out.clear();
  size_t start = 0;
  for (size_t i = 0; i <= line.size(); i++) {
    if (i == line.size() || line[i] == '\t') {
      out.push_back(line.substr(start, i - start));
      start = i + 1;
    }
  }
}

bool looks_like_header(std::string_view line, const std::vector<std::string_view>& names) {
  std::vector<std::string_view> f;
  split_tabs(line, f);
  for (auto field : f) {
    if (!field.empty() && field.back() == '\r') field.remove_suffix(1);
    for (auto name : names)
      if (field == name) return true;
  }
  return false;
}

// resolve column names against a header line.
std::vector<int> resolve_columns(std::string_view header,
                                 const std::vector<std::string_view>& names,
                                 const std::string& path) {
  std::vector<std::string_view> f;
  split_tabs(header, f);
  if (!f.empty() && !f.back().empty() && f.back().back() == '\r') f.back().remove_suffix(1);
  std::vector<int> out;
  for (auto name : names) {
    int idx = -1;
    for (int i = 0; i < (int)f.size(); i++)
      if (f[i] == name) { idx = i; break; }
    if (idx < 0)
      throw std::runtime_error(path + ": header line is missing the '" + std::string(name) +
                               "' column (every input needs a header naming its columns)");
    out.push_back(idx);
  }
  return out;
}

std::string_view nth_field(std::string_view line, int n) {
  size_t p = 0;
  for (int i = 0; i < n; i++) {
    size_t t = line.find('\t', p);
    if (t == std::string_view::npos) return {};
    p = t + 1;
  }
  size_t t = line.find('\t', p);
  return line.substr(p, t == std::string_view::npos ? std::string_view::npos : t - p);
}

std::string_view norm_chr(std::string_view chr) {
  if (chr.size() > 3 && chr.substr(0, 3) == "chr") return chr.substr(3);
  return chr;
}

bool parse_i64(std::string_view sv, int64_t& out) {
  if (sv.empty()) return false;
  int64_t val = 0;
  size_t i = (sv[0] == '-') ? 1 : 0;
  if (i >= sv.size()) return false;
  for (; i < sv.size(); i++) {
    if (sv[i] < '0' || sv[i] > '9') return false;
    val = val * 10 + (sv[i] - '0');
  }
  out = (sv[0] == '-') ? -val : val;
  return true;
}

bool normalize_motifs(std::string_view cell, std::string& out) {
  out.clear();
  size_t pos = 0;
  while (true) {
    size_t comma = cell.find(',', pos);
    std::string_view piece = cell.substr(pos, comma == std::string_view::npos ? std::string_view::npos : comma - pos);
    while (!piece.empty() && std::isspace((unsigned char)piece.front())) piece.remove_prefix(1);
    while (!piece.empty() && std::isspace((unsigned char)piece.back())) piece.remove_suffix(1);
    if (piece.empty()) return false;
    for (char c : piece)
      if (!std::strchr("ACGTacgt", c)) return false;
    if (pos) out += ',';
    out.append(piece);
    if (comma == std::string_view::npos) return true;
    pos = comma + 1;
  }
}

bool is_gzip(const std::string& path) {
  std::FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) throw std::runtime_error("cannot open: " + path);
  unsigned char magic[2] = {0, 0};
  size_t n = std::fread(magic, 1, 2, f);
  std::fclose(f);
  return n == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
}

struct RefIndex {
  bool valid = false;
  std::unordered_map<std::string, std::vector<std::pair<int64_t, int64_t>>> ranges;
};

RefIndex build_ref_index(const std::string& path, const KrefLayout& lay) {
  RefIndex ix;
  if (is_gzip(path)) return ix;

  std::FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) throw std::runtime_error("cannot open: " + path);
  std::vector<char> buf(1u << 20);
  std::string carry;
  int64_t off = 0, line_start = 0;
  std::string cur_key;
  std::vector<std::pair<int64_t, int64_t>>* cur_ranges = nullptr;

  auto note_line = [&](std::string_view line, int64_t begin, int64_t end) {
    if (begin == 0 && lay.has_header) return;
    if (line.empty()) {
      if (cur_ranges) cur_ranges->back().second = end;
      return;
    }
    std::string_view chr = norm_chr(nth_field(line, lay.chr));
    if (!cur_ranges || chr != cur_key) {
      cur_key.assign(chr);
      cur_ranges = &ix.ranges[cur_key];
      cur_ranges->push_back({begin, end});
    } else {
      cur_ranges->back().second = end;
    }
  };

  while (true) {
    size_t n = std::fread(buf.data(), 1, buf.size(), f);
    if (n == 0) break;
    size_t p = 0;
    while (p < n) {
      const char* nl = (const char*)memchr(buf.data() + p, '\n', n - p);
      if (!nl) {
        carry.append(buf.data() + p, n - p);
        break;
      }
      size_t at = (size_t)(nl - buf.data());
      int64_t line_end = off + (int64_t)at + 1;
      if (!carry.empty()) {
        carry.append(buf.data() + p, at - p);
        note_line(carry, line_start, line_end);
        carry.clear();
      } else {
        note_line(std::string_view(buf.data() + p, at - p), line_start, line_end);
      }
      p = at + 1;
      line_start = line_end;
    }
    off += (int64_t)n;
  }
  if (!carry.empty()) note_line(carry, line_start, off);  // final line without '\n'
  std::fclose(f);

  ix.valid = true;
  return ix;
}

struct TrfRow {
  std::string line;
  int64_t start, end;
};

struct ChrBlock {
  std::string chr;  // as spelled in the catalog
  std::vector<TrfRow> rows;
};

struct Catalog {
  std::string header;  // the input header line
  int n_cols = 0;
  std::vector<ChrBlock> blocks;
};

Catalog load_catalog(const std::string& path, const std::string& target, int64_t max_span) {
  GzLineReader rdr(path);
  Catalog cat;
  std::vector<ChrBlock>& blocks = cat.blocks;
  std::unordered_set<std::string> closed;
  std::string line, motifs;
  std::vector<std::string_view> f;
  uint64_t line_no = 0, n_bad_motif = 0;
  const std::string_view target_norm = norm_chr(target);

  while (rdr.next_line(line)) { line_no++; if (!line.empty()) break; }
  if (line.empty()) throw std::runtime_error(path + ": empty file");
  cat.header = line;
  const auto cols = resolve_columns(cat.header, {"chr", "tr_start", "tr_end", "tr_motif"}, path);
  const int c_chr = cols[0], c_start = cols[1], c_end = cols[2], c_motif = cols[3];
  split_tabs(cat.header, f);
  cat.n_cols = (int)f.size();
  {
    std::vector<std::string_view> suffix;
    split_tabs(OUT_SUFFIX_HEADER, suffix);
    for (auto name : suffix)
      for (auto have : f)
        if (have == name)
          throw std::runtime_error(path + ": input already has a '" + std::string(name) +
                                   "' column, which prep-catalog appends; rename it");
  }

  while (rdr.next_line(line)) {
    line_no++;
    if (line.empty()) continue;
    split_tabs(line, f);
    if ((int)f.size() != cat.n_cols)
      throw std::runtime_error(path + ": line " + std::to_string(line_no) + ": expected " +
                               std::to_string(cat.n_cols) + " columns (as in the header), got " +
                               std::to_string(f.size()));
    int64_t start, end;
    if (!parse_i64(f[c_start], start) || !parse_i64(f[c_end], end))
      throw std::runtime_error(path + ": line " + std::to_string(line_no) +
                               ": tr_start and tr_end must be integers");
    if (!normalize_motifs(f[c_motif], motifs)) { n_bad_motif++; continue; }
    if (motifs != f[c_motif]) {
      std::string fixed;
      for (size_t i = 0; i < f.size(); i++) {
        if (i) fixed += '\t';
        fixed.append(i == (size_t)c_motif ? std::string_view(motifs) : f[i]);
      }
      line = std::move(fixed);
      split_tabs(line, f);
    }
    if (end - start > max_span) continue;
    if (!target.empty() && norm_chr(f[c_chr]) != target_norm) continue;

    std::string chr(f[c_chr]);
    if (blocks.empty() || blocks.back().chr != chr) {
      if (!blocks.empty()) closed.insert(blocks.back().chr);
      if (closed.count(chr))
        throw std::runtime_error(path + ": line " + std::to_string(line_no) + ": chromosome " +
                                 chr + " appears in two separate blocks; the catalog must "
                                 "keep each chromosome contiguous");
      blocks.push_back({std::move(chr), {}});
    } else if (!blocks.back().rows.empty() && start < blocks.back().rows.back().start) {
      throw std::runtime_error(path + ": line " + std::to_string(line_no) +
                               ": rows are not sorted by start within " + chr);
    }
    blocks.back().rows.push_back({std::move(line), start, end});
  }
  if (n_bad_motif)
    std::fprintf(stderr, "[WARNING] %s: skipped %llu row(s) whose tr_motif is empty or not ACGT\n",
                 path.c_str(), (unsigned long long)n_bad_motif);
  return cat;
}

std::string read_header_line(const std::string& path) {
  GzLineReader rdr(path);
  std::string line;
  while (rdr.next_line(line))
    if (!line.empty()) return line;
  throw std::runtime_error(path + ": empty file");
}

KrefLayout read_kref_layout(const std::string& path) {
  KrefLayout L;
  const std::string first = read_header_line(path);
  const std::vector<std::string_view> names = {"chr", "position", "kmer_id", "strand"};
  if (!looks_like_header(first, names)) return L;
  const auto cols = resolve_columns(first, names, path);
  L.chr = cols[0]; L.pos = cols[1]; L.id = cols[2]; L.strand = cols[3];
  L.n_cols = 1 + std::max({L.chr, L.pos, L.id, L.strand});
  L.has_header = true;
  return L;
}

struct KRow {
  int64_t pos;
  std::string id, strand;
};

class KmerRefStream {
public:
  KmerRefStream(std::unique_ptr<LineSource> src, const std::string& path, std::string_view chr,
                const KrefLayout& lay)
      : src_(std::move(src)), path_(path), chr_norm_(norm_chr(chr)), lay_(lay) {}

  uint64_t bad_rows = 0;

  const KRow* row_before(int64_t bound) {
    if (bound <= 0) return nullptr;
    fill_past(bound - 1);
    size_t n_below = first_at_or_above(bound);
    if (n_below == 0) return nullptr;
    win_.erase(win_.begin(), win_.begin() + (n_below - 1));
    return &win_.front();
  }

  const KRow* row_after(int64_t coord) {
    fill_past(coord);
    size_t i = first_at_or_above(coord + 1);
    return i < win_.size() ? &win_[i] : nullptr;
  }

  bool has_row_between(int64_t lo, int64_t hi) {
    if (hi - KMER_LEN <= lo) return false;
    fill_past(lo);
    size_t i = first_at_or_above(lo + 1);
    return i < win_.size() && win_[i].pos < hi - KMER_LEN;
  }

private:
  std::unique_ptr<LineSource> src_;
  std::string path_;
  std::string chr_norm_;
  KrefLayout lay_;
  std::deque<KRow> win_;
  bool eof_ = false;
  int64_t last_pos_ = -1;
  uint64_t line_no_ = 0;
  std::string line_scratch_;
  std::vector<std::string_view> fields_scratch_;

  size_t first_at_or_above(int64_t x) const {
    size_t lo = 0, hi = win_.size();
    while (lo < hi) {
      size_t mid = (lo + hi) / 2;
      if (win_[mid].pos < x) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }

  void fill_past(int64_t x) {
    std::string& line = line_scratch_;
    std::vector<std::string_view>& f = fields_scratch_;
    while (!eof_ && (win_.empty() || win_.back().pos <= x)) {
      if (!src_->next_line(line)) { eof_ = true; break; }
      line_no_++;
      if (line.empty()) continue;
      split_tabs(line, f);
      if ((int)f.size() < lay_.n_cols) { bad_rows++; continue; }
      if (norm_chr(f[lay_.chr]) != chr_norm_) continue;
      int64_t pos;
      if (!parse_i64(f[lay_.pos], pos)) { bad_rows++; continue; }
      if (pos < last_pos_)
        throw std::runtime_error(path_ + ": line " + std::to_string(line_no_) +
                                 ": k-mer positions are not sorted within the chromosome");
      last_pos_ = pos;
      win_.push_back({pos, std::string(f[lay_.id]), std::string(f[lay_.strand])});
    }
  }
};

void process_chr(const ChrBlock& block, const std::string& out_header,
                 std::unique_ptr<LineSource> refs_src, const KrefLayout& kref_lay,
                 const std::string& refs_path, const std::string& piece_path,
                 std::unordered_set<std::string>& needed_ids, uint64_t& bad_rows) {
  KmerRefStream refs(std::move(refs_src), refs_path, block.chr, kref_lay);

  std::ofstream out(piece_path);
  if (!out) throw std::runtime_error("cannot write: " + piece_path);
  out << out_header << '\n';

  size_t group_begin = 0;
  int64_t group_start = 0, group_end = 0;

  auto close_group = [&](size_t group_end_idx) {
    const KRow* a = refs.row_before(group_start - KMER_LEN);
    const KRow* b = refs.row_after(group_end);
    std::string suffix;
    suffix.reserve(64);
    auto add = [&](const std::string& s) { suffix += '\t'; suffix += s; };
    add(a ? std::to_string(a->pos) : "");
    add(a ? a->id : "");
    add("NA");
    add(a ? a->strand : "");
    add(b ? std::to_string(b->pos) : "");
    add(b ? b->id : "");
    add("NA");
    add(b ? b->strand : "");
    if (a) needed_ids.insert(a->id);
    if (b) needed_ids.insert(b->id);
    for (size_t r = group_begin; r < group_end_idx; r++)
      out << block.rows[r].line << suffix << '\n';
  };

  for (size_t i = 0; i < block.rows.size(); i++) {
    const TrfRow& row = block.rows[i];
    if (i == group_begin) {
      group_start = row.start;
      group_end = row.end;
      continue;
    }
    bool merge;
    if (row.start <= group_end) {
      merge = true;
    } else {
      int64_t distance = row.start - group_end;
      merge = distance <= KMER_LEN || !refs.has_row_between(group_end, row.start);
    }
    if (merge) {
      if (row.end > group_end) group_end = row.end;
    } else {
      close_group(i);
      group_begin = i;
      group_start = row.start;
      group_end = row.end;
    }
  }
  if (group_begin < block.rows.size()) close_group(block.rows.size());

  if (!out) throw std::runtime_error("writing failed: " + piece_path);
  bad_rows = refs.bad_rows;
}

std::unordered_map<std::string, std::string>
load_kmer_seqs(const std::string& path, const std::unordered_set<std::string>& needed) {
  std::unordered_map<std::string, std::string> seqs;
  if (needed.empty()) return seqs;
  seqs.reserve(needed.size());
  GzLineReader rdr(path);
  std::string line;
  int c_id = 0, c_seq = 1;
  const std::vector<std::string_view> names = {"kmer_id", "kmer_seq"};
  bool first = true;
  std::vector<std::string_view> f;
  size_t remaining = needed.size();
  while (rdr.next_line(line)) {
    if (line.empty()) continue;
    if (first) {
      first = false;
      if (looks_like_header(line, names)) {
        const auto cols = resolve_columns(line, names, path);
        c_id = cols[0]; c_seq = cols[1];
        continue;
      }
    }
    const int n_cols = 1 + std::max(c_id, c_seq);
    split_tabs(line, f);
    if ((int)f.size() < n_cols) continue;
    std::string id(f[c_id]);
    if (!needed.count(id) || seqs.count(id)) continue;
    seqs.emplace(std::move(id), std::string(f[c_seq]));
    if (--remaining == 0) break;
  }
  return seqs;
}

void usage(std::FILE* fp, const char* prog) {
  (void)prog;
  std::fprintf(fp,
        "Usage: trumut prep-catalog -i tr_catalog.tsv -k uniqmers.tsv.gz -s uniqmer_seqs.tsv.gz -o catalog.tsv [options]\n"
        "\n"
        "Required:\n"
        "  -i, --input_catalog FILE  your tandem repeat catalog: a TSV with a header naming chr, tr_start, tr_end\n"
        "                            and tr_motif; any other columns are carried through\n"
        "  -k, --kmer_refs FILE      uniqmer positions in the reference genome (pre-built)\n"
        "  -s, --kmer_seqs FILE      uniqmer ID to sequence map (pre-built)\n"
        "  -o, --output-path FILE    output path of the uniqmer-annotated catalog, ready for trumut allele\n"
        "\n"
        "Options:\n"
        "  -t, --target STR          process specified chromosome only; 1 and chr1 both match [all]\n"
        "  -j, --jobs N              chromosomes processed in parallel [3]\n"
        "  --max-repeat-span N       drop catalog rows whose repeat spans more than N bp [15000]\n"
        "  -h, --help                show this message\n"
        "\n"
        "Documentation: https://erkim11.github.io/trumut/\n");
}

}

int prep_catalog_main(int argc, char** argv) {
  const char* catalog_path = nullptr;
  const char* seqs_path = nullptr;
  const char* refs_path = nullptr;
  const char* out_path = nullptr;
  const char* target = nullptr;
  int jobs = 3;
  int64_t max_span = 15000;

  static ko_longopt_t longopts[] = {
    { (char*)"input_catalog", ko_required_argument, 'i' },
    { (char*)"kmer_seqs",     ko_required_argument, 's' },
    { (char*)"kmer_refs",     ko_required_argument, 'k' },
    { (char*)"output-path",   ko_required_argument, 'o' },
    { (char*)"target",        ko_required_argument, 't' },
    { (char*)"jobs",          ko_required_argument, 'j' },
    { (char*)"max-repeat-span", ko_required_argument, 301 },
    { (char*)"help",          ko_no_argument,       'h' },
    { nullptr, 0, 0 }
  };
  ketopt_t k = KETOPT_INIT;
  int c;
  while ((c = ketopt(&k, argc, argv, 1, "i:s:k:o:t:j:h", longopts)) >= 0) {
    if      (c == 'i') catalog_path = k.arg;
    else if (c == 's') seqs_path    = k.arg;
    else if (c == 'k') refs_path    = k.arg;
    else if (c == 'o') out_path     = k.arg;
    else if (c == 't') target       = k.arg;
    else if (c == 'j') {
      char* end = nullptr;
      jobs = (int)std::strtol(k.arg, &end, 10);
      if (!k.arg[0] || *end || jobs < 1) { std::cerr << "[ERROR] -j must be an integer >= 1 (got '" << k.arg << "')\n"; return 2; }
    }
    else if (c == 301) {
      char* end = nullptr;
      max_span = std::strtoll(k.arg, &end, 10);
      if (!k.arg[0] || *end || max_span < 1) { std::cerr << "[ERROR] --max-repeat-span must be an integer >= 1 (got '" << k.arg << "')\n"; return 2; }
    }
    else if (c == 'h') { usage(stdout, argv[0]); return 0; }
    else { usage(stderr, argv[0]); return 2; }
  }
  if (!catalog_path || !seqs_path || !refs_path || !out_path) {
    std::cerr << "[ERROR] -i, -s, -k and -o are all required\n\n";
    usage(stderr, argv[0]);
    return 2;
  }

  try {
    Catalog cat = load_catalog(catalog_path, target ? target : "", max_span);
    std::vector<ChrBlock>& blocks = cat.blocks;
    const std::string out_header = cat.header + '\t' + OUT_SUFFIX_HEADER;
    if (blocks.empty())
      throw std::runtime_error(target
          ? std::string("chromosome '") + target + "' has no rows in " + catalog_path
          : std::string("no usable rows in ") + catalog_path);

    const fs::path out_fs(out_path);
    if (out_fs.has_parent_path()) fs::create_directories(out_fs.parent_path());
    const std::string out_dir = out_fs.has_parent_path() ? out_fs.parent_path().string() : ".";
    const std::string out_base = out_fs.filename().string();

    std::vector<std::string> pieces(blocks.size());
    for (size_t i = 0; i < blocks.size(); i++)
      pieces[i] = out_dir + "/." + out_base + "." + blocks[i].chr + ".partial";

    const KrefLayout kref_lay = read_kref_layout(refs_path);
    const RefIndex ref_index = build_ref_index(refs_path, kref_lay);
    if (!ref_index.valid && blocks.size() > 1)
      std::fprintf(stderr, "[INFO] %s is gzipped; each worker scans the whole stream (a plain copy is faster)\n", refs_path);

    const int n_jobs = (int)std::min(blocks.size(), (size_t)std::max(1, jobs));
    if (blocks.size() > 1)
      std::fprintf(stderr, "[INFO] running %zu chromosomes, %d at a time; use -j to change\n", blocks.size(), n_jobs);

    std::vector<std::unordered_set<std::string>> chr_ids(blocks.size());
    std::vector<uint64_t> chr_bad(blocks.size(), 0);
    std::vector<std::string> chr_err(blocks.size());
    std::atomic<size_t> next_chr{0};
    std::atomic<int> n_done{0};
    std::mutex log_mu;

    auto worker = [&]() {
      for (;;) {
        size_t i = next_chr.fetch_add(1);
        if (i >= blocks.size()) return;
        try {
          std::unique_ptr<LineSource> src;
          if (ref_index.valid) {
            auto it = ref_index.ranges.find(std::string(norm_chr(blocks[i].chr)));
            src = std::make_unique<RangeFileReader>(
                refs_path, it != ref_index.ranges.end()
                               ? it->second
                               : std::vector<std::pair<int64_t, int64_t>>{});
          } else {
            src = std::make_unique<GzLineReader>(refs_path);
          }
          process_chr(blocks[i], out_header, std::move(src), kref_lay, refs_path, pieces[i],
                      chr_ids[i], chr_bad[i]);
          if (blocks.size() > 1) {
            std::lock_guard<std::mutex> lg(log_mu);
            std::fprintf(stderr, "  [%2d/%zu] %s done\n", ++n_done, blocks.size(),
                         blocks[i].chr.c_str());
          }
        } catch (const std::exception& e) {
          chr_err[i] = e.what();
          std::lock_guard<std::mutex> lg(log_mu);
          std::fprintf(stderr, "  [%2d/%zu] %s FAILED: %s\n", ++n_done, blocks.size(),
                       blocks[i].chr.c_str(), e.what());
        }
      }
    };
    std::vector<std::thread> threads;
    threads.reserve(n_jobs);
    for (int i = 0; i < n_jobs; i++) threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    int n_failed = 0;
    for (const auto& e : chr_err)
      if (!e.empty()) n_failed++;
    if (n_failed) {
      for (auto& p : pieces) std::remove(p.c_str());
      throw std::runtime_error(std::to_string(n_failed) + " of " +
                               std::to_string(blocks.size()) +
                               " chromosomes failed; see the messages above");
    }

    std::unordered_set<std::string> needed_ids;
    uint64_t bad_total = 0;
    for (size_t i = 0; i < blocks.size(); i++) {
      needed_ids.insert(chr_ids[i].begin(), chr_ids[i].end());
      bad_total += chr_bad[i];
    }
    if (bad_total)
      std::fprintf(stderr, "[WARNING] skipped %llu malformed row(s) in %s\n",
                   (unsigned long long)bad_total, refs_path);

    auto seqs = load_kmer_seqs(seqs_path, needed_ids);
    if (seqs.size() < needed_ids.size())
      std::fprintf(stderr, "[WARNING] %zu of %zu kmer ids missing from %s; their sequence is NA\n",
                   needed_ids.size() - seqs.size(), needed_ids.size(), seqs_path);

    static const std::string NA = "NA";
    std::string id_scratch;
    auto seq_for = [&](std::string_view id) -> const std::string& {
      id_scratch.assign(id.begin(), id.end());
      auto it = seqs.find(id_scratch);
      return it == seqs.end() ? NA : it->second;
    };

    uint64_t rows_out = 0, rows_na = 0;
    {
      std::ofstream fout(out_path);
      if (!fout) throw std::runtime_error(std::string("cannot write: ") + out_path);
      std::string line, ob;
      ob.reserve(4096);
      std::vector<std::string_view> f;
      fout << out_header << '\n';
      for (size_t i = 0; i < blocks.size(); i++) {
        GzLineReader pin(pieces[i]);
        if (!pin.next_line(line)) continue;
        const auto cols = resolve_columns(line, {"kmerA_id", "kmerA_seq", "kmerB_id", "kmerB_seq"},
                                          pieces[i]);
        const size_t idx_a_id = cols[0], idx_a_seq = cols[1], idx_b_id = cols[2], idx_b_seq = cols[3];
        while (pin.next_line(line)) {
          split_tabs(line, f);
          const std::string& seq_a = seq_for(f[idx_a_id]);
          const std::string& seq_b = seq_for(f[idx_b_id]);
          if ((!f[idx_a_id].empty() && &seq_a == &NA) ||
              (!f[idx_b_id].empty() && &seq_b == &NA))
            rows_na++;
          ob.clear();
          for (size_t j = 0; j < f.size(); j++) {
            if (j) ob += '\t';
            if (j == idx_a_seq) ob += seq_a;
            else if (j == idx_b_seq) ob += seq_b;
            else ob.append(f[j]);
          }
          ob += '\n';
          fout.write(ob.data(), (std::streamsize)ob.size());
          rows_out++;
        }
      }
      if (!fout) throw std::runtime_error(std::string("writing failed: ") + out_path);
    }
    for (auto& p : pieces) std::remove(p.c_str());

    if (rows_na)
      std::fprintf(stderr, "[WARNING] %llu row(s) have an NA sequence; trumut allele skips them\n", (unsigned long long)rows_na);
    std::fprintf(stderr, "[INFO] wrote %llu rows to %s\n",
                 (unsigned long long)rows_out, out_path);
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "[ERROR] " << e.what() << "\n";
    return 1;
  }
}
