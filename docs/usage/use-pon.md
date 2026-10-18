# trumut use-pon

A worked example on simulated data is under [Quick start](../quick-start.md).

## Basic command

```sh
trumut use-pon [sample.tsv] [pon.tsv] [output.tsv] --input-type reads
trumut use-pon [sample.tsv] [pon.tsv] [output.tsv] --input-type assembly --ploidy [2]
```

- `sample.tsv`: sample output from `trumut allele`
- `pon.tsv`: panel of normals file from [Pre-built files](../pre-built-files.md) or `trumut build-pon`
- `output.tsv`: output path
- <code>--input-type reads&#124;assembly</code>: `reads` for sequencing reads. `assembly` for assembled contigs
- `--ploidy N`: required with `assembly`. Number of haplotype contigs per locus (e.g., `1` for a haploid assembly, `2` for diploid)

## More options

| Argument | Default | Description |
| --- | --- | --- |
| `-t chr1` | all chromosomes | process specified chromosome only, spelled as in the `trumut allele` output |
| `--clustering-thres N` | 11 | allele lengths within N bp of each other collapse into one reported allele |
| `--n-reads-thres N` | 9 | a read supports an allele if its length is within N bp of it |
| `--min-offset N` | clustering width + 1 | report an allele only if it is at least N bp from the panel range |
| `--min-supporting-reads N` | 7 | report an allele only if it has at least N supporting reads |
| `--min-motif-diff X` | 3 | report an allele only if its median motif k-mer count differs from the panel's by **more than** X |
| `--mega-events` | off | report large events only. Sets `--min-offset` 100, `--min-supporting-reads` -1 unless given explicitly (`--min-motif-diff` stays 3) |
| `--min-subclone-reads N` | 4 | call an allele `subclonal` only if another allele or density peak on the same haplotype has at least N reads |
| `--min-subclone-frac F` | 0.1 | that other allele or density peak must also hold at least this fraction of the haplotype's reads |
| `--min-reads-per-sample N` | 10; for assemblies, set by `--ploidy` | skip a locus unless the sample has at least N reads / contigs there |