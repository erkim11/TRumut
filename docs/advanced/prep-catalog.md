# trumut prep-catalog

<!-- TODO: Add paragraph describing what this does -->

## Basic command

```sh
trumut prep-catalog -i [tr_catalog.tsv] -k [uniqmers.tsv.gz] -s [uniqmer_seqs.tsv.gz] -o [catalog.tsv]
```

- `-i`: your tandem repeat catalog, a TSV with a header naming `chr`, `tr_start`, `tr_end`, and `tr_motif`. Any other columns are carried through
- `-k`: uniqmer positions in the reference genome, from [Pre-built files](../pre-built-files.md)
- `-s`: uniqmer ID to sequence map, from [Pre-built files](../pre-built-files.md)
- `-o`: output path of the uniqmer-annotated catalog, ready for `trumut allele`

## More options

| Argument | Default | Description |
| --- | --- | --- |
| `-t STR` | all | process specified chromosome only. "1" and "chr1" both match |
| `-j N` | 3 | chromosomes processed in parallel |
| `--max-repeat-span N` | 15000 | drop catalog rows whose repeat spans more than N bp |

## Output columns

| Column | Present when | Description |
| --- | --- | --- |
| *input columns* |  | every column of the input catalog, unchanged and in input order |
| `kmerA_ref_coords`, `kmerB_ref_coords` |  | 0-based reference coordinate of the two uniqmers |
| `kmerA_id`, `kmerB_id` |  | uniqmer IDs from the `-k` file |
| `kmerA_seq`, `kmerB_seq` |  | the 31-mer sequence as the reference lists it |
| `kmerA_strand`, `kmerB_strand` |  | reference strand of each uniqmer (`+`/`-`) |

A repeat group at a chromosome edge has empty coordinates, ID and strand and `NA` sequence on the missing side. `trumut allele` skips such rows.

## Example output

| chr | tr_start | tr_end | tr_motif | kmerA_ref_coords | kmerA_id | kmerA_seq | kmerA_strand | kmerB_ref_coords | kmerB_id | kmerB_seq | kmerB_strand |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 3074876 | 3074933 | CAG | 3074831 | a482910377 | GCTAAAGACAATTACATAACATACACGTCAG | + | 3074951 | a2918360054 | CACGAAACTTGTTGGCCCAGTGTGAATCGCT | - |
| chr1 | 9210440 | 9210476 | TTTA | 9210398 | a730155829 | TAAGGGTTAAGTAAGTGTGATGCATACGCCT | - | 9210492 | a1604472913 | TTACTTGCTGTGTCCACCCCATCGGACTGGC | + |
| chrX | 147912050 | 147912110 | CGG | 147912019 | x3374810259 | ATTTTTATTACACTCAGAAACAGAACTCGGG | + | 147912128 | x1029557486 | TAATTTTGACAGGTCACGCAGAGGCGCGCCC | + |