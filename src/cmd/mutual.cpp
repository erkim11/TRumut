/*
mutual.cpp -- `trcall mutual`: pairwise deviation calling across N samples.
*/

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

#include "lib/calling.h"

namespace {

struct FileGroup {
  std::unordered_map<int32_t, int32_t> freq;
  size_t n_reads = 0;
  bool covered() const { return !freq.empty(); }
  void clear() { freq.clear(); n_reads = 0; }
};

struct MutualGroup {
  std::string chr;
  int64_t min_start = std::numeric_limits<int64_t>::max();
  int64_t max_end = std::numeric_limits<int64_t>::min();
  std::string kmerA_id, kmerB_id;
  std::string loci;
  std::vector<FileGroup> per_file;

  void reset(size_t n) {
    chr.clear();
    loci.clear();
    min_start = std::numeric_limits<int64_t>::max();
    max_end = std::numeric_limits<int64_t>::min();
    per_file.assign(n, FileGroup{});
  }
};

struct Counts { int64_t a = 0, b = 0; };

struct Summary {
  std::map<std::pair<std::string, std::string>, Counts> directional;  // (test, refs_deviating)
  std::map<std::pair<std::string, std::string>, Counts> venn;         // (present, absent)
};

constexpr const char* SUMMARY_RAW_HEADER = "section\tkey1\tkey2\tcount_a\tcount_b";

void write_summary_raw(const std::string& path, const Summary& s) {
  std::ofstream out(path);
  if (!out) die("cannot write: " + path);
  out << SUMMARY_RAW_HEADER << '\n';
  for (const auto& kv : s.directional)
    out << "D\t" << kv.first.first << '\t' << kv.first.second << '\t' << kv.second.a << '\t' << kv.second.b << '\n';
  for (const auto& kv : s.venn)
    out << "V\t" << kv.first.first << '\t' << kv.first.second << '\t' << kv.second.a << '\t' << kv.second.b << '\n';
}

void merge_summary_raw(const std::string& path, Summary& s) {
  std::ifstream in(path);
  if (!in) return;
  std::string line;
  std::vector<std::string_view> f;
  if (!std::getline(in, line)) return;
  split_tsv_views(line, f);
  int c_sec = -1, c_k1 = -1, c_k2 = -1, c_a = -1, c_b = -1;
  for (int i = 0; i < (int)f.size(); ++i) {
    if      (f[i] == "section") c_sec = i;
    else if (f[i] == "key1")    c_k1 = i;
    else if (f[i] == "key2")    c_k2 = i;
    else if (f[i] == "count_a") c_a = i;
    else if (f[i] == "count_b") c_b = i;
  }
  if (c_sec < 0 || c_k1 < 0 || c_k2 < 0 || c_a < 0 || c_b < 0)
    die(path + ": summary piece header is missing a column");
  const int n_cols = 1 + std::max({c_sec, c_k1, c_k2, c_a, c_b});
  while (std::getline(in, line)) {
    split_tsv_views(line, f);
    if ((int)f.size() < n_cols) continue;
    const std::pair<std::string, std::string> key(sv_to_string(f[c_k1]), sv_to_string(f[c_k2]));
    Counts& c = (f[c_sec] == "D") ? s.directional[key] : s.venn[key];
    c.a += sv_to_i64(f[c_a]);
    c.b += sv_to_i64(f[c_b]);
  }
}

void write_summary(const std::string& path, const Summary& s) {
  std::ofstream out(path);
  if (!out) die("cannot write: " + path);
  out << "# directional: events of `test` when each input in `refs_deviating` plays the normal role\n";
  out << "section\ttest\trefs_deviating\tn_rows\tn_loci\n";
  for (const auto& kv : s.directional)
    out << "directional\t" << kv.first.first << '\t' << kv.first.second << '\t'
        << kv.second.a << '\t' << kv.second.b << '\n';
  out << "# venn: lengths carried by every input in `present` and by none in `absent`,\n"
         "# one event per locus per cluster (lengths chained within the clustering width)\n";
  out << "section\tpresent\tabsent\tn_events\tn_loci\n";
  for (const auto& kv : s.venn)
    out << "venn\t" << kv.first.first << '\t' << kv.first.second << '\t'
        << kv.second.a << '\t' << kv.second.b << '\n';
  if (!out) die("writing failed: " + path);
}

std::vector<std::string> make_labels(const std::vector<std::string>& paths) {
  std::vector<std::string> labels;
  labels.reserve(paths.size());
  for (const auto& p : paths) labels.push_back(base_name(p));
  for (size_t i = 0; i < labels.size(); ++i) {
    bool dup = false;
    for (size_t j = 0; j < labels.size(); ++j)
      if (j != i && labels[j] == labels[i]) { dup = true; break; }
    if (dup) labels[i] += "." + std::to_string(i + 1);
  }
  return labels;
}

std::string mutual_header(const std::vector<std::string>& labels) {
  std::string h =
      "chr\tregion_start\tregion_end\tkmerA_id\tkmerB_id\tloci\ttest\tallele_length\t"
      "n_reads\tn_reads_in_thres\tn_refs_deviating\trefs_deviating";
  for (const auto& l : labels) { h += "\toffset_"; h += l; }
  h += "\n";
  return h;
}

std::string join_labels(const std::vector<std::string>& labels,
                        const std::vector<int>& idx) {
  if (idx.empty()) return ".";
  std::string s;
  for (size_t i = 0; i < idx.size(); ++i) {
    if (i) s += ',';
    s += labels[idx[i]];
  }
  return s;
}

int count_length_clusters(std::vector<int32_t>& lens, int32_t threshold) {
  if (lens.empty()) return 0;
  std::sort(lens.begin(), lens.end());
  int n = 1;
  for (size_t i = 1; i < lens.size(); ++i)
    if (lens[i] - lens[i - 1] > threshold) ++n;
  return n;
}

// per-reference state of one (test, length) row.
enum class RefState : char { NotTested, Carries, Deviates };

struct RowAcc {
  std::vector<RefState> state;
  std::vector<int32_t> offset;
};

void write_group(std::ofstream& out, const MutualGroup& g,
                 const std::vector<std::string>& labels,
                 int32_t clustering_thres, int32_t n_reads_thres, int32_t min_offset,
                 size_t min_reads, Summary& summary) {
  const size_t n_files = g.per_file.size();

  std::vector<std::vector<int32_t>> lens(n_files);
  for (size_t i = 0; i < n_files; ++i) {
    const FileGroup& fg = g.per_file[i];
    if (fg.n_reads < min_reads) return;
    lens[i].reserve(fg.freq.size());
    for (const auto& kv : fg.freq) lens[i].push_back(kv.first);
  }

  std::map<std::pair<int, int32_t>, RowAcc> rows;

  for (size_t x = 0; x < n_files; ++x) {
    for (size_t y = 0; y < n_files; ++y) {
      if (y == x) continue;
      const std::vector<int32_t> dev =
          find_deviating_lengths(g.per_file[x].freq, lens[y], clustering_thres, min_offset);
      for (int32_t d : dev) {
        RowAcc& r = rows[{(int)x, d}];
        if (r.state.empty()) {
          r.state.assign(n_files, RefState::Carries);
          r.offset.assign(n_files, 0);
          r.state[x] = RefState::NotTested;
        }
        r.state[y] = RefState::Deviates;
        r.offset[y] = signed_nearest_diff_to_set(d, lens[y]);
      }
    }
  }
  if (rows.empty()) return;

  std::unordered_set<std::string> directional_seen;
  std::map<std::pair<std::string, std::string>, std::vector<int32_t>> venn_lens;

  std::string prefix = g.chr;
  prefix += '\t' + std::to_string(g.min_start);
  prefix += '\t' + std::to_string(g.max_end);
  prefix += '\t' + g.kmerA_id;
  prefix += '\t' + g.kmerB_id;
  prefix += '\t' + g.loci;

  for (const auto& kv : rows) {
    const int x = kv.first.first;
    const int32_t d = kv.first.second;
    const RowAcc& r = kv.second;

    std::vector<int> deviating, present;
    present.push_back(x);
    for (size_t z = 0; z < n_files; ++z) {
      if (r.state[z] == RefState::Deviates) deviating.push_back((int)z);
      else if (r.state[z] == RefState::Carries) present.push_back((int)z);
    }
    std::sort(present.begin(), present.end());
    const std::string dev_str = join_labels(labels, deviating);
    const std::string present_str = join_labels(labels, present);

    out << prefix
        << '\t' << labels[x]
        << '\t' << d
        << '\t' << exact_read_count(g.per_file[x].freq, d)
        << '\t' << nearby_read_count(g.per_file[x].freq, d, n_reads_thres)
        << '\t' << deviating.size()
        << '\t' << dev_str;
    for (size_t z = 0; z < n_files; ++z) {
      out << '\t';
      switch (r.state[z]) {
        case RefState::NotTested: out << "NA"; break;
        case RefState::Carries:   out << '.'; break;
        case RefState::Deviates:  out << r.offset[z]; break;
      }
    }
    out << '\n';

    {
      Counts& c = summary.directional[{labels[x], dev_str}];
      ++c.a;
      if (directional_seen.insert(labels[x] + '\t' + dev_str).second) ++c.b;
    }
    venn_lens[{present_str, dev_str}].push_back(d);
  }

  for (auto& kv : venn_lens) {
    Counts& c = summary.venn[kv.first];
    c.a += count_length_clusters(kv.second, clustering_thres);
    ++c.b;
  }
}

int run_mutual_one_chr(const std::vector<InputIndex>& ixs,
                       const std::vector<std::string>& labels,
                       const std::string& out_path, const std::string& chr,
                       int32_t clustering_thres, int32_t n_reads_thres, int32_t min_offset,
                       size_t min_reads, bool with_header,
                       const std::string& summary_raw_path) {
  std::ofstream out(out_path);
  if (!out) die("cannot write: " + out_path);
  if (with_header) out << mutual_header(labels);

  struct Stream {
    SliceReader in;
    ColIdx ci;
    std::string line;
    bool have = false;
    std::string kA, kB;
  };

  const size_t n = ixs.size();
  std::vector<Stream> st(n);
  for (size_t i = 0; i < n; ++i) {
    st[i].ci = ixs[i].ci;
    st[i].in.open(ixs[i], chr);
  }

  std::vector<std::string_view> f, lens;
  f.reserve(64);

  auto load = [&](Stream& s) {
    while (s.have) {
      if (s.line.empty()) { s.have = s.in.next(s.line); continue; }
      split_tsv_views(s.line, f);
      const int need = std::max({s.ci.chr, s.ci.region_start, s.ci.region_end, s.ci.loci,
                                 s.ci.kmerA, s.ci.kmerB, s.ci.allele_lengths});
      if ((int)f.size() <= need) { s.have = s.in.next(s.line); continue; }
      s.kA.assign(f[s.ci.kmerA].begin(), f[s.ci.kmerA].end());
      s.kB.assign(f[s.ci.kmerB].begin(), f[s.ci.kmerB].end());
      return;
    }
  };
  auto advance = [&](Stream& s) { s.have = s.in.next(s.line); load(s); };

  for (auto& s : st) { s.have = s.in.next(s.line); load(s); }

  Summary summary;
  MutualGroup g;
  while (true) {
    int lead = -1;
    for (size_t i = 0; i < n; ++i) {
      if (!st[i].have) continue;
      if (lead < 0 || cmp_key(st[i].kA, st[i].kB, st[lead].kA, st[lead].kB) < 0) lead = (int)i;
    }
    if (lead < 0) break;
    const std::string curA = st[lead].kA, curB = st[lead].kB;

    g.reset(n);
    g.kmerA_id = curA;
    g.kmerB_id = curB;

    for (size_t i = 0; i < n; ++i) {
      Stream& s = st[i];
      while (s.have && s.kA == curA && s.kB == curB) {
        split_tsv_views(s.line, f);
        if (g.chr.empty())  g.chr  = sv_to_string(f[s.ci.chr]);
        if (g.loci.empty()) g.loci = sv_to_string(f[s.ci.loci]);
        const int64_t start = sv_to_i64(f[s.ci.region_start]);
        const int64_t end   = sv_to_i64(f[s.ci.region_end]);
        if (start < g.min_start) g.min_start = start;
        if (end   > g.max_end)   g.max_end   = end;
        // one row is one k-mer pair
        split_list(f[s.ci.allele_lengths], lens);
        for (auto sv : lens) {
          const int32_t L = (int32_t)sv_to_i64(sv);
          ++g.per_file[i].freq[L];
          ++g.per_file[i].n_reads;
        }
        advance(s);
      }
    }

    if (g.chr.empty() || g.min_start == std::numeric_limits<int64_t>::max() ||
        g.max_end == std::numeric_limits<int64_t>::min())
      continue;

    write_group(out, g, labels, clustering_thres, n_reads_thres, min_offset, min_reads, summary);
  }

  if (!out) die("writing failed: " + out_path);
  if (!summary_raw_path.empty()) write_summary_raw(summary_raw_path, summary);
  return 0;
}

int run_mutual(const std::vector<std::string>& paths, const std::string& out_path,
               const std::string& target_chr, int32_t fixed_clustering_thres,
               int32_t fixed_n_reads_thres, int32_t fixed_min_offset,
               size_t min_reads, const std::string& summary_path, int n_jobs_max) {
  if (out_path.empty()) die("output path is empty");
  const std::string outdir = parent_dir(out_path);
  ensure_dir_exists(outdir);

  const std::vector<std::string> labels = make_labels(paths);
  const int32_t clustering_thres = fixed_clustering_thres >= 0 ? fixed_clustering_thres : 3;
  const int32_t n_reads_thres = fixed_n_reads_thres >= 0 ? fixed_n_reads_thres : clustering_thres;
  const int32_t min_offset = fixed_min_offset >= 0 ? fixed_min_offset : clustering_thres + 1;
  const bool want_summary = !summary_path.empty();
  const std::string out_base = base_name(out_path);

  std::vector<InputIndex> ixs(paths.size());
  for (size_t i = 0; i < paths.size(); ++i)
    build_input_index(paths[i], "mutual" + std::to_string(i + 1), outdir, ixs[i]);
  auto cleanup_temps = [&]() {
    for (const auto& ix : ixs) if (ix.is_temp) std::remove(ix.path.c_str());
  };

  std::vector<std::string> chrs;
  if (!target_chr.empty()) {
    chrs.push_back(target_chr);
  } else {
    for (const auto& ix : ixs)
      for (const auto& r : ix.ranges)
        if (std::find(chrs.begin(), chrs.end(), r.chr) == chrs.end()) chrs.push_back(r.chr);
  }
  if (chrs.empty()) { cleanup_temps(); die("no data rows in any input"); }

  auto finish_summary = [&](const std::vector<std::string>& raw_paths) {
    if (!want_summary) return;
    Summary merged;
    for (const auto& p : raw_paths) { merge_summary_raw(p, merged); std::remove(p.c_str()); }
    write_summary(summary_path, merged);
  };

  if (chrs.size() == 1) {
    const std::string raw = want_summary ? outdir + "/." + out_base + ".summary.partial" : std::string();
    const int rc = run_mutual_one_chr(ixs, labels, out_path, chrs[0], clustering_thres, n_reads_thres, min_offset,
                                      min_reads, /*with_header=*/true, raw);
    cleanup_temps();
    if (rc == 0) finish_summary(raw.empty() ? std::vector<std::string>{} : std::vector<std::string>{raw});
    return rc;
  }

  std::vector<std::string> pieces(chrs.size()), raws(chrs.size());
  for (size_t i = 0; i < chrs.size(); ++i) {
    pieces[i] = outdir + "/." + out_base + "." + chrs[i] + ".partial";
    if (want_summary) raws[i] = pieces[i] + ".summary";
  }

  const int n_jobs = (int)std::min(chrs.size(), (size_t)std::max(1, n_jobs_max));
  std::fprintf(stderr, "[INFO] running %zu chromosomes, %d at a time; use -j to change\n", chrs.size(), n_jobs);

  size_t next = 0;
  int running = 0, n_failed = 0, n_done = 0;
  std::vector<pid_t> pids(chrs.size(), -1);
  std::vector<char> ok(chrs.size(), 0);

  while (next < chrs.size() || running > 0) {
    while (running < n_jobs && next < chrs.size()) {
      std::fflush(nullptr);
      const pid_t pid = fork();
      if (pid < 0) {
        std::fprintf(stderr, "[ERROR] cannot fork a worker for %s\n", chrs[next].c_str());
        ++n_failed; ++next;
        continue;
      }
      if (pid == 0)
        _exit(run_mutual_one_chr(ixs, labels, pieces[next], chrs[next], clustering_thres, n_reads_thres, min_offset,
                                 min_reads, /*with_header=*/false, raws[next]));
      pids[next] = pid;
      ++running; ++next;
    }

    int status = 0;
    const pid_t done_pid = waitpid(-1, &status, 0);
    if (done_pid < 0) {
      if (errno == EINTR) continue;
      break;
    }
    for (size_t i = 0; i < pids.size(); ++i) {
      if (pids[i] != done_pid) continue;
      pids[i] = -1;
      --running; ++n_done;
      const int rc = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
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
    fout << mutual_header(labels);
    for (size_t i = 0; i < chrs.size(); ++i) {
      if (ok[i]) {
        std::ifstream pin(pieces[i], std::ios::binary);
        if (pin && pin.peek() != std::char_traits<char>::eof()) fout << pin.rdbuf();
      }
      std::remove(pieces[i].c_str());
    }
    if (!fout) die("writing failed: " + out_path);
  }

  if (want_summary) {
    std::vector<std::string> good;
    for (size_t i = 0; i < chrs.size(); ++i) {
      if (ok[i]) good.push_back(raws[i]);
      else std::remove(raws[i].c_str());
    }
    finish_summary(good);
  }

  if (n_failed)
    std::fprintf(stderr, "[ERROR] %zu of %zu chromosomes failed; %s is incomplete (%zu succeeded), see messages above\n",
                 (size_t)n_failed, chrs.size(), out_path.c_str(), chrs.size() - (size_t)n_failed);
  return n_failed ? 1 : 0;
}

void mutual_usage(std::ostream& os, const char* prog) {
  (void)prog;
  os <<
    "Usage: trcall mutual output.tsv -i sample1.tsv -i sample2.tsv [-i sample3.tsv ...] [options]\n"
    "\n"
    "Positional:\n"
    "  output.tsv       output path\n"
    "\n"
    "Required:\n"
    "  -i FILE                   trcall allele output; at least two must be given\n"
    "\n"
    "Options:\n"
    "  -t chr1                   process specified chromosome only, spelled as in the trcall allele output [all]\n"
    "  -j N                      chromosomes processed in parallel [3]\n"
    "  --clustering-thres N      allele lengths within N bp of each other collapse into one reported allele [3]\n"
    "  --n-reads-thres N         a read supports an allele if its length is within N bp of it [= clustering width]\n"
    "  --min-offset N            report a length only if it is at least N bp from every other input's lengths [clustering width + 1]\n"
    "  --min-reads-per-sample N  skip a locus unless each sample has at least N reads there [10]\n"
    "  --summary PATH            write a two-section summary: row counts per test sample and reference set\n"
    "                            (directional), and allele counts per present/absent sample combination (venn)\n"
    "  -h, --help                show this message\n"
    "\n"
    "Documentation: https://erkim11.github.io/trcall/\n";
}

}  // namespace

int mutual_main(int argc, char** argv) {
  std::vector<std::string> inputs, positionals;
  std::string target_chr, summary_path;
  int n_jobs = 3;   // -j: chromosomes at once
  int32_t fixed_clustering_thres = -1;
  int32_t fixed_n_reads_thres = -1;
  int32_t fixed_min_offset = -1;
  size_t min_reads = MIN_UNIQUE_READS_PER_INPUT;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      inputs.push_back(argv[++i]);
    } else if (arg == "-t" && i + 1 < argc) {
      target_chr = argv[++i];
    } else if (arg == "-j" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], n_jobs)) return 2;
      if (n_jobs < 1) { std::cerr << "[ERROR] -j must be at least 1\n"; return 2; }
    } else if ((arg == "--clustering-thres" || arg == "--n-reads-thres") && i + 1 < argc) {
      int v;
      if (!parse_int_arg(arg.c_str(), argv[++i], v)) return 2;
      if (v < 0) {
        std::cerr << "[ERROR] " << arg << " must be a non-negative integer\n";
        return 2;
      }
      if (arg == "--clustering-thres") fixed_clustering_thres = v;
      else                             fixed_n_reads_thres = v;
    } else if (arg == "--min-offset" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], fixed_min_offset)) return 2;
      if (fixed_min_offset < 0) {
        std::cerr << "[ERROR] --min-offset must be a non-negative integer\n";
        return 2;
      }
    } else if (arg == "--min-reads-per-sample" && i + 1 < argc) {
      int mr;
      if (!parse_int_arg(arg.c_str(), argv[++i], mr)) return 2;
      if (mr < 0) {
        std::cerr << "[ERROR] --min-reads-per-sample must be a non-negative integer\n";
        return 2;
      }
      min_reads = (size_t)mr;
    } else if (arg == "--summary" && i + 1 < argc) {
      summary_path = argv[++i];
    } else if (arg == "-h" || arg == "--help") {
      mutual_usage(std::cout, argv[0]);
      return 0;
    } else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "[ERROR] unknown option: " << arg << "\n";
      mutual_usage(std::cerr, argv[0]);
      return 2;
    } else {
      positionals.push_back(arg);
    }
  }

  if (positionals.empty()) {
    mutual_usage(std::cerr, argv[0]);
    return 2;
  }
  const std::string out_path = positionals.front();
  inputs.insert(inputs.end(), positionals.begin() + 1, positionals.end());

  if (inputs.size() < 2) {
    std::cerr << "[ERROR] at least two inputs are required\n";
    mutual_usage(std::cerr, argv[0]);
    return 2;
  }

  return run_mutual(inputs, out_path, target_chr, fixed_clustering_thres, fixed_n_reads_thres,
                    fixed_min_offset, min_reads, summary_path, n_jobs);
}
