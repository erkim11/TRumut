# trumut build-pon

<!-- TODO: Add paragraph describing what this does -->

## Basic command

```sh
trumut build-pon [alleles_dir/] [pon.tsv] --input-type reads
trumut build-pon [alleles_dir/] [pon.tsv] --input-type assembly --ploidy 2
```

- `[alleles_dir/]`: directory holding one `trumut allele` output per normal sample, named `*.alleles.tsv`
- `[pon.tsv]`: output path of the panel of normals
- <code>--input-type reads&#124;assembly</code>: `reads` for sequencing reads. `assembly` for assembled contigs
- `--ploidy N`: required with `assembly`. Number of haplotype contigs per locus (e.g., `1` for a haploid assembly, `2` for diploid)

### With digests

Digesting each allele output first is optional. A digest is much smaller in file size and merges to an identical panel, so it pays off when the allele outputs live on slow storage or the panel is rebuilt often.

```sh
trumut build-pon --digest sample1.alleles.tsv sample1.pon-digest.tsv   # one sample
trumut build-pon --digest alleles_dir/   # every *.alleles.tsv in the directory

# merge as before
trumut build-pon alleles_dir/ pon.tsv --input-type reads
```

- `sample1.pon-digest.tsv`: one digest per sample, written next to the allele output or into an optional output directory given as a third argument
- the merge uses a sample's digest when one exists and falls back to the raw allele output otherwise

## More options

| Argument | Default | Description |
| --- | --- | --- |
| `--pair-percent-support F` | 0.98 | fraction of ALL samples that must cover a pair for it to be kept |
| `--distance-percent-support F` | 0.0975 | fraction of the COVERING samples that must carry a length for it to count as supported |
| `--min-reads-per-sample N` | 10; for assemblies, set by `--ploidy` | a sample "covers" a pair only with at least N reads/contigs there |
| `--suffix S` | .alleles.tsv | input filename suffix |
| `--jobs J` | 3 | merge tasks at once, or files digested at once in `--digest` directory mode |
| `--split-factor N` | by chromosome size | cut every chromosome into exactly N kmerA_id sub-slices for parallelization |
| `--digest` | off | reduces raw outputs to smaller digests that merge to the same panel |

## Output columns

**Panel table**

| Column | Present when | Description |
| --- | --- | --- |
| `chr`, `kmerA_id`, `kmerB_id` |  | copied from the `trumut allele` output |
| `range_min`, `range_max` |  | smallest / largest allele length carried by at least `--distance-percent-support` of the covering samples |
| `range_min_n_motif_kmers`, `range_max_n_motif_kmers` |  | median across samples of each sample's median number of motif k-mers\* at that length |

**Digest** (`--digest`)

| Column | Present when | Description |
| --- | --- | --- |
| `chr`, `kmerA_id`, `kmerB_id` |  | copied from the `trumut allele` output |
| `n_reads` |  | number of distinct reads at the pair in that sample |
| `allele_lengths` |  | distinct allele lengths (list); unlike the allele output, one item per distinct length, not per read |
| `motif_kmers` |  | median number of motif k-mers\* per allele length (list) |

<p class="table-note">* An approximate measure of how much of the allele is made of the repeat motif.</p>

Columns marked (list) hold one comma-separated item per read, in the same order across columns.

## Example output

**Panel table**

| chr | kmerA_id | kmerB_id | range_min | range_max | range_min_n_motif_kmers | range_max_n_motif_kmers |
| --- | --- | --- | --- | --- | --- | --- |
| chr1 | a482910377 | a2918360054 | 80 | 107 | 45 | 72 |
| chr1 | a730155829 | a1604472913 | 55 | 75 | 24 | 43.5 |
| chr9 | a2207731158 | a689043327 | NA | NA | NA | NA |

**Digest**

| chr | kmerA_id | kmerB_id | n_reads | allele_lengths | motif_kmers |
| --- | --- | --- | --- | --- | --- |
| chr1 | a482910377 | a2918360054 | 31 | 89,92,87,...,90 | 51,54,49,...,52 |
| chr1 | a730155829 | a1604472913 | 36 | 63,67,65,...,71 | 30,34,32,...,38 |
| chr9 | a2207731158 | a689043327 | 4 | 412,419,416,422 | 287,294,291,297 |