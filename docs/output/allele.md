# trumut allele

## Output columns

| Column | Present when | Description |
| --- | --- | --- |
| `chr` |  | chromosome of the uniqmer pair |
| `region_start` |  | leftmost repeat start coordinate between the two uniqmers (0-based) |
| `region_end` |  | rightmost repeat end coordinate between the two uniqmers |
| `loci` |  | catalog tandem repeats between the two uniqmers, as `start-end:motif` |
| `query_names` | `--save-query-names` | read name, or contig name for assemblies (list) |
| `kmerA_id`, `kmerB_id` |  | the two uniqmer IDs |
| `kmerA_ref_coords`, `kmerB_ref_coords` | `--save-kmer-coords` | reference start of each uniqmer (0-based) |
| `n_reads` |  | number of reads (or contigs) in which both uniqmers were found |
| `allele_lengths` |  | number of bases between the two uniqmers (list) |
| `n_motif_kmers` |  | number of k-mers in the allele sequence that match the repeat motif\* (list) |
| `allele_seqs` | `--save-allele-seq` | sequence between the two uniqmers (list) |
| `HP` |  | haplotype tag, `1` or `2` (list); `.` if untagged |
| `PS` |  | phase set (list); `.` whenever `HP` is `.` |

<p class="table-note">* An approximate measure of how much of the allele is made of the repeat motif.</p>

Columns with an empty "Present when" cell are always written. Columns marked (list) hold one comma-separated item per read, in the same order across columns. `.` marks a missing value.

## Example output

| chr | region_start | region_end | loci | query_names | kmerA_id | kmerB_id | kmerA_ref_coords | kmerB_ref_coords | n_reads | allele_lengths | n_motif_kmers | allele_seqs | HP | PS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 3074876 | 3074933 | 3074876-3074933:CAG | m84011/1201/ccs,m84011/4417/ccs,m84011/9020/ccs,...,m84011/9977/ccs | a482910377 | a2918360054 | 3074831 | 3074951 | 31 | 89,89,92,...,89 | 51,51,54,...,51 | GATTCGCTAGGTACCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTCCGATAGCTTGACGGAT,GATTCGCTAGGTACCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTCCGATAGCTTGACGGAT,GATTCGCTAGGTACCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTCCGATAGCTTGACGGAT,...,GATTCGCTAGGTACCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTCCGATAGCTTGACGGAT | 1,2,1,...,2 | 3050112,3050112,3050112,...,3050112 |
| chr1 | 9210440 | 9210476 | 9210440-9210476:TTTA | m84011/0302/ccs,m84011/5560/ccs,m84011/7781/ccs,...,m84011/8890/ccs | a730155829 | a1604472913 | 9210398 | 9210492 | 36 | 63,67,63,...,67 | 30,34,30,...,34 | CCTGAAGTCGATTTATTTATTTATTTATTTATTTATTTATTTATTTAAGGCTTACGATCGTTC,CCTGAAGTCGATTTATTTATTTATTTATTTATTTATTTATTTATTTATTTAAGGCTTACGATCGTTC,CCTGAAGTCGATTTATTTATTTATTTATTTATTTATTTATTTATTTAAGGCTTACGATCGTTC,...,CCTGAAGTCGATTTATTTATTTATTTATTTATTTATTTATTTATTTATTTAAGGCTTACGATCGTTC | 1,.,2,...,2 | 9188004,.,9188004,...,9188004 |
| chr2 | 1188020 | 1188110 | 1188020-1188061:AT,1188070-1188110:GT | m84011/2210/ccs,m84011/3388/ccs,m84011/6104/ccs,...,m84011/7015/ccs | a1049382267 | a377258106 | 1187970 | 1188122 | 26 | 121,121,121,...,121 | 68,68,68,...,68 | TGCAGTCAAGGCTTACGATATATATATATATATATATATATATATATATATATATATATCGGTCAGCGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCGTAGGCTTACA,TGCAGTCAAGGCTTACGATATATATATATATATATATATATATATATATATATATATATCGGTCAGCGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCGTAGGCTTACA,TGCAGTCAAGGCTTACGATATATATATATATATATATATATATATATATATATATATATCGGTCAGCGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCGTAGGCTTACA,...,TGCAGTCAAGGCTTACGATATATATATATATATATATATATATATATATATATATATATCGGTCAGCGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCGTAGGCTTACA | .,.,.,...,. | .,.,.,...,. |

## Optional output files

Three optional files record what did not make it into the main output. Each is written only when its flag is given.

### `-w, --wrong-kmer-hits`

Reads / contigs in which both uniqmers were found once, but on mismatched strands or in the wrong order, so the length between them is not trusted. Same columns as the main output.

| chr | region_start | region_end | loci | kmerA_id | kmerB_id | n_reads | allele_lengths | n_motif_kmers | HP | PS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 3074876 | 3074933 | 3074876-3074933:CAG | a482910377 | a2918360054 | 2 | 89,89 | 51,51 | 1,. | 3050112,. |
| chr5 | 146878727 | 146878759 | 146878727-146878759:CTG | a1520497731 | a44062815 | 1 | 71 | 26 | 2 | 146850990 |
| chr11 | 2884420 | 2884480 | 2884420-2884480:GGA | a2231705644 | a985170328 | 3 | 78,78,81 | 54,54,57 | .,.,. | .,.,. |

### `-m, --many-kmer-hits`

Uniqmer pairs that were found more than once in a single read / contig, so no single length could be measured. The columns depend on the input type.

**Reads input**: one row per read and pair.

| Column | Description |
| --- | --- |
| `query_name` | read name |
| `chr`, `region_start`, `region_end` | the pair's region, from the catalog |
| `kmerA_id`, `kmerB_id` | the two uniqmer IDs |
| `kmerA_hit_count`, `kmerB_hit_count` | times each uniqmer was found in the read |

| query_name | chr | region_start | region_end | kmerA_id | kmerB_id | kmerA_hit_count | kmerB_hit_count |
| --- | --- | --- | --- | --- | --- | --- | --- |
| m84011/2731/ccs | chr1 | 3074876 | 3074933 | a482910377 | a2918360054 | 2 | 1 |
| m84011/5108/ccs | chr8 | 118371245 | 118371301 | a1806613252 | a3149920067 | 1 | 3 |
| m84011/6620/ccs | chr8 | 118371245 | 118371301 | a1806613252 | a3149920067 | 2 | 2 |

**Assembly input**: one row per pair, judged over the whole assembly.

| Column | Description |
| --- | --- |
| `chr`, `region_start`, `region_end` | the pair's region, from the catalog |
| `kmerA_id`, `kmerB_id` | the two uniqmer IDs |
| `kmerA_contig`, `kmerB_contig` | the first contig each uniqmer was found on |
| `kmerA_location`, `kmerB_location` | position of that hit on that contig (0-based) |
| `kmerA_hit_count`, `kmerB_hit_count` | total hits of each uniqmer across all contigs |

| chr | region_start | region_end | kmerA_id | kmerB_id | kmerA_contig | kmerA_location | kmerA_hit_count | kmerB_contig | kmerB_location | kmerB_hit_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 3074876 | 3074933 | a482910377 | a2918360054 | h1tg000004l | 3061207 | 3 | h1tg000004l | 3061327 | 3 |
| chr8 | 118371245 | 118371301 | a1806613252 | a3149920067 | h1tg000031l | 9945112 | 1 | h2tg000018l | 4420571 | 4 |
| chr17 | 44907520 | 44907571 | a3308811109 | a1274459938 | h2tg000012l | 21883540 | 2 | h2tg000012l | 21883640 | 1 |

### `-d, --debug-file`

One row per uniqmer pair in the catalog, saying whether it reached the main output and, if not, why.

| Column | Description |
| --- | --- |
| `kmerA_id`, `kmerB_id` | the two uniqmer IDs |
| `dropped` | `yes` if the pair is absent from the main output, else `no` |
| `reason` | why the pair was dropped, as a short phrase: one or both uniqmers not found, found too many times, never found together, or found on mismatched strands or in the wrong order. Empty when the pair was not dropped |

| kmerA_id | kmerB_id | dropped | reason |
| --- | --- | --- | --- |
| a482910377 | a2918360054 | no |  |
| a1806613252 | a3149920067 | yes | kmer(s) hit more than once in a read |
| a3308811109 | a1274459938 | yes | kmers never found in the same read |
