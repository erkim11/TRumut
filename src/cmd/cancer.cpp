/*
cancer.cpp -- `trumut cancer`: reports alleles in a tumor that are absent from its matched normal
*/

#include <iostream>
#include <string>
#include "lib/calling.h"

static void cancer_usage(std::ostream& os, const char* prog) {
  (void)prog;
  os <<
    "Usage: trumut cancer tumor.tsv normal.tsv output.tsv --input-type reads|assembly [options]\n"
    "\n"
    "Positional:\n"
    "  tumor.tsv        tumor output from trumut allele\n"
    "  normal.tsv       matched normal output from trumut allele\n"
    "  output.tsv       output path\n"
    "\n"
    "Required:\n"
    "  --input-type STR          reads for sequencing reads, assembly for assembled contigs\n"
    "  --ploidy N                required with assembly: haplotype contigs per locus (1 haploid, 2 diploid)\n"
    "\n"
    "Options:\n"
    "  -t chr1                   process specified chromosome only, spelled as in the trumut allele output [all]\n"
    "  -j N                      chromosomes processed in parallel [3]\n"
    "  --clustering-thres N      allele lengths within N bp of each other collapse into one reported allele [11]\n"
    "  --n-reads-thres N         a read supports an allele if its length is within N bp of it [9]\n"
    "  --min-offset N            report a tumor allele only if it is at least N bp from every normal allele [clustering width + 1]\n"
    "  --min-supporting-reads N  report a tumor allele only if it has at least N supporting reads [7]\n"
    "  --min-motif-diff X        report a tumor allele only if its median motif k-mer count differs from the normal's by more than X [3]\n"
    "  --mega-events             report large events only: sets --min-offset 100, --min-supporting-reads -1 unless given explicitly (--min-motif-diff stays 3)\n"
    "  --min-subclone-reads N    call an allele subclonal only if another allele or density peak on the same haplotype has at least N reads [4]\n"
    "  --min-subclone-frac F     that other allele or density peak must also hold at least this fraction of the haplotype's reads [0.1]\n"
    "  --min-tagged-frac F       use haplotype-aware event types only if at least this fraction of the locus's reads carry an HP tag [0.5]\n"
    "  --min-tagged-reads N      use haplotype-aware event types only if at least N HP-tagged reads support the allele [2]\n"
    "  --min-gap-kmers N         count a non-motif gap toward n_nonmotif_kmers_in only if it holds at least N k-mers [10]\n"
    "  --min-gap-frac F          count a non-motif gap only if at least this fraction of its k-mers exceed the reference [0.95]\n"
    "  --min-reads-per-sample N  skip a locus unless tumor and normal each have at least N reads / contigs there [10; for assemblies, set by --ploidy]\n"
    "  -h, --help                show this message\n"
    "\n"
    "Documentation: https://erkim11.github.io/trumut/\n";
}

int cancer_main(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      cancer_usage(std::cout, argv[0]);
      return 0;
    }
  }
  if (argc < 4) {
    cancer_usage(std::cerr, argv[0]);
    return 2;
  }
  const std::string tumor_path  = argv[1];
  const std::string normal_path = argv[2];
  const std::string out_path    = argv[3];

  std::string target_chr;        // empty means all chromosomes, in parallel
  int32_t fixed_clustering_thres = -1;  // -1 means use the default (11)
  int32_t fixed_n_reads_thres = -1;     // -1 means use the default (9)
  int32_t fixed_min_offset = -1;        // -1 means use the clustering width
  int input_type = 0;   // --input-type: 1 reads, -1 assembly; required
  int32_t fixed_min_reads = -1;  // -1 means the --input-type default (10 or 2)
  int ploidy = 0;                // --ploidy: 0 means not given (required for assemblies)
  int n_jobs = 3;                // -j: chromosomes at once
  RowFilters filters;   // --min-supporting-reads / --min-motif-diff
  bool saw_supporting = false, saw_motif_diff = false;
  bool mega_events = false;      // --mega-events: shifts the two defaults above

  // parse optional arguments
  for (int i = 4; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-t" && i + 1 < argc) {
      target_chr = argv[++i];
    } else if (arg == "-j" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], n_jobs)) return 2;
      if (n_jobs < 1) {
        std::cerr << "[ERROR] -j must be at least 1\n";
        return 2;
      }
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
    } else if (arg == "--min-supporting-reads" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], filters.min_supporting_reads)) return 2;
      saw_supporting = true;
    } else if (arg == "--min-motif-diff" && i + 1 < argc) {
      if (!parse_double_arg(arg.c_str(), argv[++i], filters.min_motif_diff)) return 2;    
      saw_motif_diff = true;
    } else if (arg == "--min-subclone-reads" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], filters.min_subclone_reads)) return 2;
      if (filters.min_subclone_reads < 0) {
        std::cerr << "[ERROR] --min-subclone-reads must be a non-negative integer\n";
        return 2;
      }
    } else if (arg == "--min-subclone-frac" && i + 1 < argc) {
      if (!parse_double_arg(arg.c_str(), argv[++i], filters.min_subclone_frac)) return 2;
      if (filters.min_subclone_frac < 0.0 || filters.min_subclone_frac > 1.0) {
        std::cerr << "[ERROR] --min-subclone-frac must be between 0 and 1\n";
        return 2;
      }
    } else if (arg == "--min-tagged-frac" && i + 1 < argc) {
      if (!parse_double_arg(arg.c_str(), argv[++i], filters.min_tagged_frac)) return 2;
      if (filters.min_tagged_frac < 0.0 || filters.min_tagged_frac > 1.0) {
        std::cerr << "[ERROR] --min-tagged-frac must be between 0 and 1\n";
        return 2;
      }
    } else if (arg == "--min-tagged-reads" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], filters.min_tagged_reads)) return 2;
      if (filters.min_tagged_reads < 0) {
        std::cerr << "[ERROR] --min-tagged-reads must be a non-negative integer\n";
        return 2;
      }
    } else if (arg == "--min-gap-kmers" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], filters.min_gap_kmers)) return 2;
      if (filters.min_gap_kmers < 0) {
        std::cerr << "[ERROR] --min-gap-kmers must be a non-negative integer\n";
        return 2;
      }
    } else if (arg == "--min-gap-frac" && i + 1 < argc) {
      if (!parse_double_arg(arg.c_str(), argv[++i], filters.min_gap_frac)) return 2;
      if (filters.min_gap_frac < 0.0 || filters.min_gap_frac > 1.0) {
        std::cerr << "[ERROR] --min-gap-frac must be between 0 and 1\n";
        return 2;
      }
    } else if (arg == "--mega-events") {
      mega_events = true;
    } else if (arg == "--input-type" && i + 1 < argc) {
      input_type = parse_input_type(argv[++i]);
      if (input_type == 0) {
        std::cerr << "[ERROR] --input-type must be 'reads' or 'assembly' (got '" << argv[i] << "')\n";
        return 2;
      }
    } else if (arg == "--min-reads-per-sample" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], fixed_min_reads)) return 2;
      if (fixed_min_reads < 0) {
        std::cerr << "[ERROR] --min-reads-per-sample must be a non-negative integer\n";
        return 2;
      }
    } else if (arg == "--ploidy" && i + 1 < argc) {
      if (!parse_int_arg(arg.c_str(), argv[++i], ploidy)) return 2;
      if (ploidy < 1) {
        std::cerr << "[ERROR] --ploidy must be a positive integer (1 haploid, 2 diploid, ...)\n";
        return 2;
      }
    } else {
      std::cerr << "[ERROR] unknown option: " << arg << "\n";
      cancer_usage(std::cerr, argv[0]);
      return 2;
    }
  }

  if (input_type == 0) {
    std::cerr << "[ERROR] " << INPUT_TYPE_REQUIRED_MSG;
    return 2;
  }
  if (mega_events) {   // large events only
    if (fixed_min_offset < 0) fixed_min_offset = 100;
    if (!saw_supporting)      filters.min_supporting_reads = -1;   // disabled
    if (!saw_motif_diff)      filters.min_motif_diff = 3.0;
  }
  const bool assembly = input_type == -1;
  const size_t min_reads = resolve_min_reads(assembly, ploidy, fixed_min_reads, "[ERROR] ");
  if (min_reads == MIN_READS_INVALID) return 2;

  return run_deviation_calling(tumor_path, {normal_path}, out_path, target_chr,
                               fixed_clustering_thres, fixed_n_reads_thres, fixed_min_offset,
                               {"tumor", "normal"}, filters, min_reads, assembly, n_jobs);
}
