# trumut cancer

A worked example on simulated data is under [Quick start](../quick-start.md).

## Basic command

```sh
trumut cancer [tumor.tsv] [normal.tsv] [output.tsv] --input-type reads
trumut cancer [tumor.tsv] [normal.tsv] [output.tsv] --input-type assembly --ploidy [2]
```

- `[tumor.tsv]`: tumor output from `trumut allele`
- `[normal.tsv]`: matched normal output from `trumut allele`
- `[output.tsv]`: output path
- <code>--input-type reads&#124;assembly</code>: `reads` for sequencing reads. `assembly` for assembled contigs
- `--ploidy N`: required with `assembly`. Number of haplotype contigs per locus (e.g., `1` for a haploid assembly, `2` for diploid)

## More options

| Argument | Default | Description |
| --- | --- | --- |
| `-t chr1` | all chromosomes | process specified chromosome only, spelled as in the `trumut allele` output |
| `-j N` | 3 | chromosomes processed in parallel |
| `--clustering-thres N` | 11 | allele lengths within N bp of each other collapse into one reported allele |
| `--n-reads-thres N` | 9 | a read supports an allele if its length is within N bp of it |
| `--min-offset N` | clustering width + 1 | report a tumor allele only if it is at least N bp from every normal allele |
| `--min-supporting-reads N` | 7 | report a tumor allele only if it has at least N supporting reads |
| `--min-motif-diff X` | 3 | report a tumor allele only if its median motif k-mer count differs from the normal's by **more than** X |
| `--mega-events` | off | report large events only. Sets `--min-offset` 100, `--min-supporting-reads` -1 unless given explicitly (`--min-motif-diff` stays 3) |
| `--min-reads-per-sample N` | 10; for assemblies, set by `--ploidy` | skip a locus unless tumor and normal each have at least N reads / contigs there |
| `--min-subclone-reads N` | 4 | call a tumor allele `subclonal` only if another allele or density peak on the same haplotype has at least N reads |
| `--min-subclone-frac F` | 0.1 | that other allele or density peak must also hold at least this fraction of the haplotype's reads |
| `--min-gap-kmers N` | 10 | count a non-motif insertion toward `n_nonmotif_kmers_in` column value only if it holds at least N k-mers |
| `--min-gap-frac F` | 0.95 | count a non-motif gap only if at least this fraction of its k-mers are newly inserted |