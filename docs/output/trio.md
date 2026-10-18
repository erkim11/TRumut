# trumut trio

## Output columns

| Column | Present when | Description |
| --- | --- | --- |
| `chr`, `region_start`, `region_end`, `kmerA_id`, `kmerB_id`, `loci` |  | copied from the `trumut allele` output |
| `allele_length` |  | child allele length that is absent from the parents |
| `offset` |  | child `allele_length` minus the nearest parent allele length |
| `n_motif_kmers_ref_median` |  | median number of motif k-mers\* over parent reads at that nearest length |
| `n_motif_kmers_test_median` |  | median number of motif k-mers\* over child reads at `allele_length` |
| `n_nonmotif_kmers_in` | `--save-allele-seq` in both inputs | number of non-motif k-mers (i.e., non-repetitive sequences) inserted within the repeat |
| `n_reads` | reads input only | number of child reads at exactly `allele_length` |
| `n_reads_in_thres` | reads input only | number of child reads supporting the `allele_length` (within `--n-reads-thres` of `allele_length` and at least `--min-offset` from every parent length) |
| `event_type` |  | clonality of the mutation event. `clonal`, `subclonal`, or `private_phased` if the locus is sufficiently haplotagged (`--min-tagged-frac`, `--min-tagged-reads`). Otherwise `private_unphased` or `not_private_unphased`. Always `clonal` for assemblies |
| `phase_counts` | `HP` column in the test input | `HP=1:HP=2:untagged` counts among reads that support `allele_length` |

<p class="table-note">* An approximate measure of how much of the allele is made of the repeat motif.</p>

Columns with an empty "Present when" cell are always written.

!!! warning "event_type without haplotags"
    Without enough haplotagged reads (see `--min-tagged-frac` and `--min-tagged-reads`), `clonal` and `subclonal` cannot be told apart and are both reported as `not_private_unphased`. `private_unphased` is likewise less reliable than `private_phased`.
    
## Example output

| chr | region_start | region_end | kmerA_id | kmerB_id | loci | allele_length | offset | n_motif_kmers_ref_median | n_motif_kmers_test_median | n_nonmotif_kmers_in | n_reads | n_reads_in_thres | event_type | phase_counts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 3074876 | 3074933 | a482910377 | a2918360054 | 3074876-3074933:CAG | 39 | -18 | 51 | 33 | 0 | 16 | 19 | clonal | 19:0:0 |
| chr12 | 6936716 | 6936775 | a1682044917 | a3105628741 | 6936716-6936775:CAG | 98 | 40 | 53 | 92 | 0 | 2 | 7 | private_unphased | 0:0:7 |
| chr3 | 63912685 | 63912790 | a2588104376 | a941730268 | 63912685-63912740:CAG,63912752-63912790:CAA | 141 | 36 | 94 | 130 | 5 | 11 | 13 | subclonal | 13:0:0 |
