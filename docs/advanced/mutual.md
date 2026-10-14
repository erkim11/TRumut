# trcall mutual

<!-- TODO: Add paragraph describing what this does -->

## Basic command

```sh
trcall mutual [output.tsv] -i [sample1.tsv] -i [sample2.tsv] ...
```

- `[output.tsv]`: output path
- `-i`: `trcall allele` output. At least two must be given

## More options

| Argument | Default | Description |
| --- | --- | --- |
| `-t chr1` | all | process specified chromosome only, spelled as in the `trcall allele` output |
| `-j N` | 3 | chromosomes processed in parallel |
| `--clustering-thres N` | 3 | allele lengths within N bp of each other collapse into one reported allele |
| `--n-reads-thres N` | = clustering width | a read supports an allele if its length is within N bp of it |
| `--min-offset N` | clustering width + 1 | report a length only if it is at least N bp from every other input's lengths |
| `--min-reads-per-sample N` | 10 | skip a locus unless each sample has at least N reads there |
| `--summary path` | none | write a two-section summary: row counts per `test` sample and reference set (directional), and allele counts per present/absent sample combination (venn) |

## Output columns

| Column | Present when | Description |
| --- | --- | --- |
| `chr`, `region_start`, `region_end`, `kmerA_id`, `kmerB_id`, `loci` |  | copied from the `trcall allele` output |
| `test` |  | label (file basename) of the input the allele was found in |
| `allele_length` |  | allele length in the `test` sample |
| `n_reads` |  | number of `test` sample reads at exactly `allele_length` |
| `n_reads_in_thres` |  | number of `test` sample reads supporting the `allele_length` (within `--n-reads-thres` of `allele_length`) |
| `n_refs_deviating` |  | number of other samples that lack the allele |
| `refs_deviating` |  | comma-separated labels of the samples that lack the allele |
| `offset_<sample>` (one column per input) |  | `allele_length` minus that sample's nearest allele length. `.` if the sample carries the allele. `NA` in the `test` sample's own column |

## Example output

| chr | region_start | region_end | kmerA_id | kmerB_id | loci | test | allele_length | n_reads | n_reads_in_thres | n_refs_deviating | refs_deviating | offset_s1.tsv | offset_s2.tsv | offset_s3.tsv |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 3074876 | 3074933 | a482910377 | a2918360054 | 3074876-3074933:CAG | s1.tsv | 107 | 14 | 16 | 2 | s2.tsv,s3.tsv | NA | 18 | 21 |
| chr1 | 9210440 | 9210476 | a730155829 | a1604472913 | 9210440-9210476:TTTA | s2.tsv | 55 | 9 | 11 | 1 | s3.tsv | . | NA | -8 |
| chr2 | 1188020 | 1188110 | a1049382267 | a377258106 | 1188020-1188061:AT,1188070-1188110:GT | s3.tsv | 133 | 12 | 12 | 2 | s1.tsv,s2.tsv | 12 | 12 | NA |

**Summary file**

**Directional** counts output rows, grouped by the `test` sample and the set of samples it deviates from. The row with `test = s1.tsv` and `refs_deviating = s3.tsv` is exactly what `trcall cancer s1.tsv s3.tsv` would report. Because every row has one `test` sample, an allele shared by s1 and s2 but absent from s3 appears twice: once with s1 as `test` and once with s2 as `test`.

| section | test | refs_deviating | n_rows | n_loci |
| --- | --- | --- | --- | --- |
| directional | s1.tsv | s3.tsv | 158 | 150 |
| directional | s2.tsv | s3.tsv | 154 | 147 |
| directional | s1.tsv | s2.tsv,s3.tsv | 265 | 251 |

**Venn** counts alleles, grouped by which samples carry them (`present`) and which lack them (`absent`). It merges the output rows to count deviating events detected per sample combination. Each row is one region of a Venn diagram of the samples.

| section | present | absent | n_events | n_loci |
| --- | --- | --- | --- | --- |
| venn | s1.tsv,s2.tsv | s3.tsv | 146 | 139 |
| venn | s1.tsv | s2.tsv,s3.tsv | 251 | 240 |
| venn | s2.tsv,s3.tsv | s1.tsv | 172 | 163 |