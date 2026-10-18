/*
use-pon.cpp -- `trumut use-pon`: reports alleles outside the range seen in a panel of normals
*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "lib/calling.h"
#include "third_party/ketopt.h"

struct PonEntry {
  int32_t min_sup = 0, max_sup = 0;
  bool sup_valid = false;
  double med_motif_min = std::numeric_limits<double>::quiet_NaN();
  double med_motif_max = std::numeric_limits<double>::quiet_NaN();
};
using PonMap = std::unordered_map<std::string, PonEntry>;

static inline std::string make_pair_key(std::string_view a, std::string_view b) {
  std::string key;
  key.reserve(a.size() + 1 + b.size());
  key.append(a);
  key += '\t';
  key.append(b);
  return key;
}

static double parse_double_or_nan(std::string_view sv) {
  if (sv.empty()) return std::numeric_limits<double>::quiet_NaN();
  try {
    return std::stod(std::string(sv));
  } catch (...) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

static std::string format_double(double v) {
  char buf[32];
  int n = std::snprintf(buf, sizeof buf, "%g", v);
  return std::string(buf, (size_t)n);
}

static bool parse_i32(std::string_view sv, int32_t* out) {
  if (sv.empty()) return false;
  int64_t v = 0;
  size_t i = (sv[0] == '-') ? 1 : 0;
  if (i >= sv.size()) return false;
  for (; i < sv.size(); ++i) {
    if (sv[i] < '0' || sv[i] > '9') return false;
    v = v * 10 + (sv[i] - '0');
    if (v > std::numeric_limits<int32_t>::max()) return false;
  }
  *out = (int32_t)(sv[0] == '-' ? -v : v);
  return true;
}

// load the trumut build-pon table.
static PonMap load_pon(const char* path) {
  LineReader in;   // plain or gzipped
  if (!in.open(path)) die(std::string("cannot open PON file: ") + path);

  PonMap pon;
  std::string line;

  if (!in.getline(line)) die("PON file is empty");

  std::vector<std::string_view> fields;
  fields.reserve(8);
  split_tsv_views(line, fields);

  int col_kmerA = -1, col_kmerB = -1, col_min_sup = -1, col_max_sup = -1;
  int col_med_min = -1, col_med_max = -1;
  for (int i = 0; i < (int)fields.size(); ++i) {
    if      (fields[i] == "kmerA_id")                   col_kmerA   = i;
    else if (fields[i] == "kmerB_id")                   col_kmerB   = i;
    else if (fields[i] == "range_min")                  col_min_sup = i;
    else if (fields[i] == "range_max")                  col_max_sup = i;
    else if (fields[i] == "range_min_n_motif_kmers")  col_med_min = i;
    else if (fields[i] == "range_max_n_motif_kmers")  col_med_max = i;
  }
  if (col_kmerA   < 0) die("PON missing column: kmerA_id");
  if (col_kmerB   < 0) die("PON missing column: kmerB_id");
  if (col_min_sup < 0) die("PON missing column: range_min (is this a trumut build-pon output?)");
  if (col_max_sup < 0) die("PON missing column: range_max");
  if (col_med_min < 0) die("PON missing column: range_min_n_motif_kmers");
  if (col_med_max < 0) die("PON missing column: range_max_n_motif_kmers");

  const int min_cols = std::max({col_kmerA, col_kmerB, col_min_sup, col_max_sup,
                                 col_med_min, col_med_max}) + 1;

  while (in.getline(line)) {
    if (line.empty()) continue;
    split_tsv_views(line, fields);
    if ((int)fields.size() < min_cols) continue;

    PonEntry entry;
    entry.sup_valid = parse_i32(fields[col_min_sup], &entry.min_sup) &&
                      parse_i32(fields[col_max_sup], &entry.max_sup);
    entry.med_motif_min = parse_double_or_nan(fields[col_med_min]);
    entry.med_motif_max = parse_double_or_nan(fields[col_med_max]);

    pon[make_pair_key(fields[col_kmerA], fields[col_kmerB])] = entry;
  }

  return pon;
}

struct PonLookup {
  bool in_pon = false;      // pair present in the PoN at all
  bool sup_valid = false;
  bool survive = false;
  int32_t offset = 0;
  double normal_motif = std::numeric_limits<double>::quiet_NaN();
};

static PonLookup lookup_pon(const PonMap& pon, const std::string& key,
                            int32_t dist, int32_t min_offset) {
  PonLookup r;
  auto pair_it = pon.find(key);
  if (pair_it == pon.end()) return r;
  r.in_pon = true;

  const PonEntry& e = pair_it->second;
  if (!e.sup_valid) return r;
  r.sup_valid = true;

  r.survive = dist <= e.min_sup - min_offset || dist >= e.max_sup + min_offset;

  const int32_t off_min = dist - e.min_sup;
  const int32_t off_max = dist - e.max_sup;
  if (std::abs(off_min) <= std::abs(off_max)) {
    r.offset = off_min;
    r.normal_motif = e.med_motif_min;
  } else {
    r.offset = off_max;
    r.normal_motif = e.med_motif_max;
  }
  return r;
}

struct MosaicGroupData {
  std::string chr;
  int64_t min_start = std::numeric_limits<int64_t>::max();
  int64_t max_end   = std::numeric_limits<int64_t>::min();
  std::string kmerA_id, kmerB_id, loci;

  std::unordered_map<int32_t,int32_t> surv_freq;
  std::unordered_map<int32_t,PonLookup> pon_by_dist;

  std::unordered_map<int32_t,int32_t> test_freq;
  std::unordered_map<int32_t,HpCounts> test_hp_freq;
  std::unordered_map<int32_t, std::vector<double>> motif_by_len;

  size_t n_test_reads = 0;
};

// write one row per deviating (clustered) distance in the group.
static void write_mosaic_rows(
    std::ostream& out,
    const MosaicGroupData& g,
    const std::vector<int32_t>& deviating,
    int32_t n_reads_thres,
    const RowFilters& filters,
    bool with_phase,
    bool assembly)
{
  const LengthPredicate is_deviating =
      [&](int32_t L) { return g.pon_by_dist.count(L) > 0; };
  for (int32_t d : deviating) {
    auto pl_it = g.pon_by_dist.find(d);
    const PonLookup* pl = (pl_it != g.pon_by_dist.end()) ? &pl_it->second : nullptr;
    const int32_t offset_val = pl ? pl->offset : 0;

    std::string normal_motif = "nan";
    if (pl && !std::isnan(pl->normal_motif)) normal_motif = format_double(pl->normal_motif);

    const int32_t num_exact  = exact_read_count(g.test_freq, d);
    const int32_t num_nearby = nearby_read_count(g.test_freq, d, n_reads_thres, is_deviating);
    const std::string phase = with_phase
        ? phase_counts_for_length(g.test_hp_freq, d, n_reads_thres) : std::string();
    const char* event_type = assembly
        ? "clonal"
        : event_type_for_length(g.test_hp_freq, d, n_reads_thres, is_deviating, filters);

    double tumor_motif_v = std::numeric_limits<double>::quiet_NaN();
    if (auto it = g.motif_by_len.find(d); it != g.motif_by_len.end()) {
      std::vector<double> vals = it->second;
      tumor_motif_v = median_of(vals);
    }
    const double normal_motif_v = pl ? pl->normal_motif : std::numeric_limits<double>::quiet_NaN();
    if (!filters.keep(assembly, num_nearby, tumor_motif_v, normal_motif_v)) continue;
    const std::string tumor_motif = format_median(tumor_motif_v, "nan");

    out << g.chr << '\t' << g.min_start << '\t' << g.max_end << '\t'
        << g.kmerA_id << '\t' << g.kmerB_id << '\t' << g.loci
        << '\t' << d << '\t' << offset_val << '\t'
        << normal_motif << '\t' << tumor_motif;
    if (!assembly) out << '\t' << num_exact << '\t' << num_nearby;
    out << '\t' << event_type;
    if (with_phase) out << '\t' << phase;
    out << '\n';
  }
}

static std::string mosaic_header(bool with_phase, bool assembly) {
  std::string h =
      "chr\tregion_start\tregion_end\tkmerA_id\tkmerB_id\tloci\tallele_length\toffset\t"
      "n_motif_kmers_ref_median\tn_motif_kmers_test_median";
  if (!assembly) h += "\tn_reads\tn_reads_in_thres";
  h += "\tevent_type";
  if (with_phase) h += "\tphase_counts";
  h += "\n";
  return h;
}

static void call_mosaic(const char* tumor_path, const PonMap& pon,
                        const char* output_path, const char* target_chr,
                        int32_t clustering_thres, int32_t n_reads_thres,
                        int32_t min_offset,
                        size_t min_reads, const RowFilters& filters, bool assembly) {
  LineReader in;   // plain or gzipped (trumut allele -o *.gz)
  if (!in.open(tumor_path)) die(std::string("cannot open tumor file: ") + tumor_path);

  const bool to_stdout = !output_path || std::string_view(output_path) == "-";
  std::ofstream file_out;
  if (!to_stdout) {
    file_out.open(output_path);
    if (!file_out) die(std::string("cannot open output file: ") + output_path);
  }
  std::ostream& out = to_stdout ? std::cout : static_cast<std::ostream&>(file_out);

  std::string header;
  if (!in.getline(header)) die("tumor file is empty");

  ColIdx ci;
  if (!parse_header_indices(header, ci)) die("missing required columns in header: " + std::string(tumor_path));
  const bool with_phase = ci.hp >= 0;
  const int max_col = std::max({ci.chr, ci.region_start, ci.region_end, ci.loci, ci.kmerA, ci.kmerB,
                                ci.allele_lengths, ci.n_motif_kmers});

  out << mosaic_header(with_phase, assembly);

  std::string line;
  std::vector<std::string_view> f;
  f.reserve(64);
  RowReads rr;
  bool have_line = in.getline(line);

  while (have_line) {
    split_tsv_views(line, f);
    if ((int)f.size() <= max_col) { have_line = in.getline(line); continue; }
    if (target_chr && f[ci.chr] != std::string_view(target_chr)) { have_line = in.getline(line); continue; }

    const std::string curA = sv_to_string(f[ci.kmerA]);
    const std::string curB = sv_to_string(f[ci.kmerB]);

    MosaicGroupData g;
    g.kmerA_id = curA;
    g.kmerB_id = curB;

    while (have_line) {
      split_tsv_views(line, f);
      if ((int)f.size() <= max_col) { have_line = in.getline(line); continue; }

      std::string a = sv_to_string(f[ci.kmerA]);
      std::string b = sv_to_string(f[ci.kmerB]);
      if (a != curA || b != curB) break;

      if (g.chr.empty())  g.chr  = sv_to_string(f[ci.chr]);
      if (g.loci.empty()) g.loci = sv_to_string(f[ci.loci]);
      int64_t st = sv_to_i64(f[ci.region_start]);
      int64_t en = sv_to_i64(f[ci.region_end]);
      if (st < g.min_start) g.min_start = st;
      if (en > g.max_end)   g.max_end   = en;

      const std::string key = make_pair_key(curA, curB);
      expand_row_reads(f, ci, rr);
      for (size_t i = 0; i < rr.n; ++i) {
        const int32_t dist = (int32_t)sv_to_i64(rr.lens[i]);
        ++g.n_test_reads;
        const PonLookup pl = lookup_pon(pon, key, dist, min_offset);
        const bool survive = pl.in_pon && pl.sup_valid && pl.survive;

        ++g.test_freq[dist];
        g.motif_by_len[dist].push_back(parse_motif_count(RowReads::at(rr.motifs, i)));
        const std::string_view h = RowReads::at(rr.hp, i);
        if      (h == "1") ++g.test_hp_freq[dist].hp1;
        else if (h == "2") ++g.test_hp_freq[dist].hp2;
        else               ++g.test_hp_freq[dist].na;

        if (survive) {
          ++g.surv_freq[dist];
          g.pon_by_dist[dist] = pl;
        }
      }

      have_line = in.getline(line);
    }

    if (g.n_test_reads < min_reads) continue; // too few supporting reads for this loci

    if (!g.surv_freq.empty() && g.min_start != std::numeric_limits<int64_t>::max()) {
      std::vector<LenCand> cands;
      cands.reserve(g.surv_freq.size());
      for (const auto& kv : g.surv_freq) cands.push_back({kv.first, kv.second});

      std::vector<int32_t> deviating = greedy_cluster_lengths(std::move(cands), clustering_thres);
      if (!deviating.empty())
        write_mosaic_rows(out, g, deviating, n_reads_thres, filters, with_phase, assembly);
    }
  }
}

static void use_pon_usage(std::ostream& os) {
  os <<
    "Usage: trumut use-pon sample.tsv pon.tsv output.tsv --input-type reads|assembly [options]\n"
    "\n"
    "Positional:\n"
    "  sample.tsv       sample output from trumut allele\n"
    "  pon.tsv          panel of normals (pre-built, or from trumut build-pon)\n"
    "  output.tsv       output path\n"
    "\n"
    "Required:\n"
    "  --input-type STR          reads for sequencing reads, assembly for assembled contigs\n"
    "  --ploidy N                required with assembly: haplotype contigs per locus (1 haploid, 2 diploid)\n"
    "\n"
    "Options:\n"
    "  -t, --target STR          process specified chromosome only, spelled as in the trumut allele output [all]\n"
    "  --clustering-thres N      allele lengths within N bp of each other collapse into one reported allele [7]\n"
    "  --n-reads-thres N         a read supports an allele if its length is within N bp of it [9]\n"
    "  --min-offset N            report an allele only if it is at least N bp from the panel range [clustering width + 1]\n"
    "  --min-supporting-reads N  report an allele only if it has at least N supporting reads [7]\n"
    "  --min-motif-diff X        report an allele only if its median motif k-mer count differs from the panel's by more than X [3]\n"
    "  --mega-events             report large events only: sets --min-offset 100, --min-supporting-reads -1 unless given explicitly (--min-motif-diff stays 3)\n"
    "  --min-subclone-reads N    call an allele subclonal only if another allele or density peak on the same haplotype has at least N reads [4]\n"
    "  --min-subclone-frac F     that other allele or density peak must also hold at least this fraction of the haplotype's reads [0.1]\n"
    "  --min-tagged-frac F       use haplotype-aware event types only if at least this fraction of the locus's reads carry an HP tag [0.5]\n"
    "  --min-tagged-reads N      use haplotype-aware event types only if at least N HP-tagged reads support the allele [2]\n"
    "  --min-reads-per-sample N  skip a locus unless the sample has at least N reads / contigs there [10; for assemblies, set by --ploidy]\n"
    "  -h, --help                show this message\n"
    "\n"
    "Documentation: https://erkim11.github.io/trumut/\n";
}

int use_pon_main(int argc, char* argv[]) {
  const char* tumor_path  = nullptr;
  const char* pon_path    = nullptr;
  const char* output_path = nullptr;
  const char* target_chr  = nullptr;
  int32_t clustering_thres = 7;
  int32_t n_reads_thres = -1;  // -1: default to 9
  int32_t min_offset = -1;     // -1: default to clustering_thres
  size_t min_reads = 10;
  int fixed_min_reads = -1;   // -m: -1 means not given
  int ploidy = 0;             // --ploidy: 0 means not given
  RowFilters filters;   // --min-supporting-reads / --min-motif-diff
  bool saw_supporting = false, saw_motif_diff = false;
  bool mega_events = false;   // --mega-events: shifts two defaults
  int input_type = 0;   // --input-type: 1 reads, -1 assembly; required

  static ko_longopt_t longopts[] = {
    { (char*)"target",           ko_required_argument, 't' },
    { (char*)"clustering-thres", ko_required_argument, 305 },
    { (char*)"n-reads-thres",    ko_required_argument, 306 },
    { (char*)"min-reads-per-sample", ko_required_argument, 309 },
    { (char*)"min-offset",       ko_required_argument, 302 },
    { (char*)"min-supporting-reads", ko_required_argument, 303 },
    { (char*)"min-motif-diff",   ko_required_argument, 304 },
    { (char*)"input-type",       ko_required_argument, 301 },
    { (char*)"ploidy",           ko_required_argument, 307 },
    { (char*)"mega-events",      ko_no_argument,       308 },
    { (char*)"min-subclone-reads", ko_required_argument, 310 },
    { (char*)"min-subclone-frac", ko_required_argument, 311 },
    { (char*)"min-tagged-frac", ko_required_argument, 312 },
    { (char*)"min-tagged-reads", ko_required_argument, 313 },
    { (char*)"help",             ko_no_argument,       'h' },
    { 0, 0, 0 }
  };

  ketopt_t k = KETOPT_INIT;
  int c;
  while ((c = ketopt(&k, argc, argv, 1, "t:h", longopts)) >= 0) {
    if      (c == 't') target_chr  = k.arg;
    else if (c == 305) {
      if (!parse_int_arg("--clustering-thres", k.arg, clustering_thres)) return 1;
      if (clustering_thres < 0) { std::cerr << "[ERROR] --clustering-thres must be a non-negative integer\n"; return 1; }
    }
    else if (c == 306) {
      if (!parse_int_arg("--n-reads-thres", k.arg, n_reads_thres)) return 1;
      if (n_reads_thres < 0) { std::cerr << "[ERROR] --n-reads-thres must be a non-negative integer\n"; return 1; }
    }
    else if (c == 309) {
      int mr;
      if (!parse_int_arg("--min-reads-per-sample", k.arg, mr)) return 1;
      if (mr < 0) { std::cerr << "[ERROR] --min-reads-per-sample must be a non-negative integer\n"; return 1; }
      fixed_min_reads = mr;
    }
    else if (c == 307) {
      if (!parse_int_arg("--ploidy", k.arg, ploidy)) return 1;
      if (ploidy < 1) { std::cerr << "[ERROR] --ploidy must be a positive integer (1 haploid, 2 diploid, ...)\n"; return 1; }
    }
    else if (c == 302) {
      if (!parse_int_arg("--min-offset", k.arg, min_offset)) return 1;
      if (min_offset < 0) { std::cerr << "[ERROR] --min-offset must be a non-negative integer\n"; return 1; }
    }
    else if (c == 301) {
      input_type = parse_input_type(k.arg);
      if (input_type == 0) { std::cerr << "[ERROR] --input-type must be 'reads' or 'assembly' (got '" << k.arg << "')\n"; return 1; }
    }
    else if (c == 303) { if (!parse_int_arg("--min-supporting-reads", k.arg, filters.min_supporting_reads)) return 1; saw_supporting = true; }
    else if (c == 304) { if (!parse_double_arg("--min-motif-diff", k.arg, filters.min_motif_diff)) return 1; saw_motif_diff = true; }
    else if (c == 308) mega_events = true;
    else if (c == 310) {
      if (!parse_int_arg("--min-subclone-reads", k.arg, filters.min_subclone_reads)) return 1;
      if (filters.min_subclone_reads < 0) { std::cerr << "[ERROR] --min-subclone-reads must be a non-negative integer\n"; return 1; }
    }
    else if (c == 311) {
      if (!parse_double_arg("--min-subclone-frac", k.arg, filters.min_subclone_frac)) return 1;
      if (filters.min_subclone_frac < 0.0 || filters.min_subclone_frac > 1.0) { std::cerr << "[ERROR] --min-subclone-frac must be between 0 and 1\n"; return 1; }
    }
    else if (c == 312) {
      if (!parse_double_arg("--min-tagged-frac", k.arg, filters.min_tagged_frac)) return 1;
      if (filters.min_tagged_frac < 0.0 || filters.min_tagged_frac > 1.0) { std::cerr << "[ERROR] --min-tagged-frac must be between 0 and 1\n"; return 1; }
    }
    else if (c == 313) {
      if (!parse_int_arg("--min-tagged-reads", k.arg, filters.min_tagged_reads)) return 1;
      if (filters.min_tagged_reads < 0) { std::cerr << "[ERROR] --min-tagged-reads must be a non-negative integer\n"; return 1; }
    }
    else if (c == 'h') { use_pon_usage(std::cout); return 0; }
    else if (c == '?') { std::cerr << "[ERROR] unknown option\n"; use_pon_usage(std::cerr); return 1; }
    else if (c == ':') { std::cerr << "[ERROR] missing argument\n"; use_pon_usage(std::cerr); return 1; }
  }

  if (argc - k.ind != 3) {
    if (argc - k.ind < 3)
      std::cerr << "[ERROR] use-pon takes three files: tumor.tsv pon.tsv output.tsv\n";
    else
      std::cerr << "[ERROR] unexpected extra argument '" << argv[k.ind + 3]
                << "'; use-pon takes three files: tumor.tsv pon.tsv output.tsv\n";
    use_pon_usage(std::cerr);
    return 1;
  }
  tumor_path  = argv[k.ind];
  pon_path    = argv[k.ind + 1];
  output_path = argv[k.ind + 2];

  if (input_type == 0) {
    std::cerr << "[ERROR] " << INPUT_TYPE_REQUIRED_MSG;
    return 1;
  }
  const bool assembly = input_type == -1;
  min_reads = resolve_min_reads(assembly, ploidy, fixed_min_reads, "[ERROR] ");
  if (min_reads == MIN_READS_INVALID) return 1;

  if (mega_events) {   // large events only
    if (min_offset < 0)  min_offset = 100;
    if (!saw_supporting) filters.min_supporting_reads = -1;   // disabled
    if (!saw_motif_diff) filters.min_motif_diff = 3.0;
  }
  if (n_reads_thres < 0) n_reads_thres = 9;
  if (min_offset < 0) min_offset = clustering_thres + 1;

  PonMap pon = load_pon(pon_path);
  call_mosaic(tumor_path, pon, output_path, target_chr, clustering_thres, n_reads_thres, min_offset, min_reads,
              filters, assembly);

  return 0;
}
