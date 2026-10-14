# trcall cancer

## Output columns

| Column | Present when | Description |
| --- | --- | --- |
| `chr`, `region_start`, `region_end`, `kmerA_id`, `kmerB_id`, `loci` |  | copied from the `trcall allele` output |
| `allele_length` |  | tumor allele length that is absent from the matched normal |
| `offset` |  | tumor `allele_length` minus the nearest normal allele length |
| `n_motif_kmers_ref_median` |  | median number of motif k-mers\* over normal reads at that nearest length |
| `n_motif_kmers_test_median` |  | median number of motif k-mers\* over tumor reads at `allele_length` |
| `n_nonmotif_kmers_in` | `--save-allele-seq` in both inputs | number of non-motif k-mers (i.e., non-repetitive sequences) inserted within the repeat |
| `n_reads` | reads input only | number of tumor reads at exactly `allele_length` |
| `n_reads_in_thres` | reads input only | number of tumor reads supporting the `allele_length` (within `--n-reads-thres` of `allele_length`) |
| `event_type` |  | clonality of the mutation event. `clonal`, `subclonal`, or `private_phased` if the window has haplotagged reads. Otherwise `private_unphased` or `not_private_unphased`. Always `clonal` for assemblies |
| `phase_counts` | `HP` column in the test input | `HP=1:HP=2:untagged` counts among reads that support `allele_length` |

<p class="table-note">* An approximate measure of how much of the allele is made of the repeat motif.</p>

Columns with an empty "Present when" cell are always written.

!!! warning "event_type without haplotags"
    Without haplotagged reads, `clonal` and `subclonal` cannot be told apart and are both reported as `not_private_unphased`. `private_unphased` is likewise less reliable than `private_phased`.

## Example output

| chr | region_start | region_end | kmerA_id | kmerB_id | loci | allele_length | offset | n_motif_kmers_ref_median | n_motif_kmers_test_median | n_nonmotif_kmers_in | n_reads | n_reads_in_thres | event_type | phase_counts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr4 | 3076604 | 3076660 | a915027744 | a2263791458 | 3076604-3076660:CAG | 33 | -24 | 51 | 27 | 0 | 12 | 12 | clonal | 12:0:0 |
| chr9 | 27573482 | 27573546 | a2207731158 | a689043327 | 27573482-27573546:GGGGCC | 340 | 276 | 58 | 330 | 18 | 7 | 10 | private_phased | 1:9:0 |
| chr16 | 87637880 | 87637999 | a1937426053 | a608215794 | 87637880-87637919:AAAG,87637941-87637999:AAAAG | 168 | 49 | 105 | 152 | 11 | 8 | 10 | clonal | 0:10:0 |
