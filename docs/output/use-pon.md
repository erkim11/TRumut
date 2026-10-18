# trumut use-pon

## Output columns

| Column | Present when | Description |
| --- | --- | --- |
| `chr`, `region_start`, `region_end`, `kmerA_id`, `kmerB_id`, `loci` |  | copied from the `trumut allele` output |
| `allele_length` |  | sample allele length that is not within the panel range |
| `offset` |  | sample `allele_length` minus the closer panel bound (`range_min` or `range_max`); negative below the range, positive above |
| `n_motif_kmers_ref_median` |  | the panel's `range_min_n_motif_kmers` or `range_max_n_motif_kmers` at that closer bound; `nan` when the panel value is `NA` |
| `n_motif_kmers_test_median` |  | median number of motif k-mers\* over sample reads at `allele_length` |
| `n_reads` | reads input only | number of sample reads at exactly `allele_length` |
| `n_reads_in_thres` | reads input only | number of sample reads supporting the `allele_length` (within `--n-reads-thres` of `allele_length`) |
| `event_type` |  | clonality of the mutation event. `clonal`, `subclonal`, or `private_phased` if the locus is sufficiently haplotagged (`--min-tagged-frac`, `--min-tagged-reads`). Otherwise `private_unphased` or `not_private_unphased`. Always `clonal` for assemblies |
| `phase_counts` | `HP` column in the test input | `HP=1:HP=2:untagged` counts among reads that support `allele_length` |

<p class="table-note">* An approximate measure of how much of the allele is made of the repeat motif.</p>

Columns with an empty "Present when" cell are always written.

!!! warning "event_type without haplotags"
    Without enough haplotagged reads (see `--min-tagged-frac` and `--min-tagged-reads`), `clonal` and `subclonal` cannot be told apart and are both reported as `not_private_unphased`. `private_unphased` is likewise less reliable than `private_phased`.
    
## Example output

| chr | region_start | region_end | kmerA_id | kmerB_id | loci | allele_length | offset | n_motif_kmers_ref_median | n_motif_kmers_test_median | n_reads | n_reads_in_thres | event_type | phase_counts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 3074876 | 3074933 | a482910377 | a2918360054 | 3074876-3074933:CAG | 93 | 30 | 51 | 87 | 8 | 9 | subclonal | 2:7:0 |
| chr7 | 55019017 | 55019062 | a1355219086 | a2496073815 | 55019017-55019062:TTTA | 12 | -20 | 27 | 8 | 6 | 6 | not_private_unphased | 0:0:6 |
| chr22 | 46191235 | 46191330 | a3061557129 | a827904416 | 46191235-46191270:ATTCT,46191284-46191330:ATTCC | 136 | 41 | 84 | 121 | 7 | 9 | clonal | 4:5:0 |
