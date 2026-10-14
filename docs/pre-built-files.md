# Pre-built files

All files below are in the following Zenodo record: [https://doi.org/10.5281/zenodo.22903449](https://doi.org/10.5281/zenodo.22903449).

## Downloading

Download from the record page in a browser or from the command line with `wget` or `curl` using the file links on that page:

```sh
wget https://zenodo.org/records/22903449/files/hg38.trf-b7-b7-p500.uniqmers.tsv
wget https://zenodo.org/records/22903449/files/HPRC.pon.09012026.tsv
wget https://zenodo.org/records/22903449/files/hs38.ukm579.s02.sorted.txt.gz
wget https://zenodo.org/records/22903449/files/ukm579.s02.txt.gz
```

## Uniqmer-annotated tandem repeat catalog

**hg38.trf-b7-b7-p500.uniqmers.tsv**

The catalog `trcall allele` reads, made by `trcall prep-catalog` from a genome-wide tandem repeat catalog.

| Column | Description |
| --- | --- |
| `chr` | chromosome |
| `tr_start` | repeat start (0-based) |
| `tr_end` | repeat end |
| `tr_motif` | repeat motif |
| `kmerA_ref_coords`, `kmerB_ref_coords` | 0-based reference coordinate of the two uniqmers |
| `kmerA_id`, `kmerB_id` | uniqmer IDs |
| `kmerA_seq`, `kmerB_seq` | the 31-mer sequence as the reference lists it  |
| `kmerA_strand`, `kmerB_strand` | reference strand of each uniqmer (`+`/`-`) |

## Panel of normals

**HPRC.pon.09012026.tsv**

The panel `trcall use-pon` compares against, made by `trcall build-pon` from normal samples. One row per uniqmer pair.

| Column | Description |
| --- | --- |
| `chr` | chromosome |
| `kmerA_id`, `kmerB_id` | uniqmer IDs |
| `range_min`, `range_max` | smallest / largest allele length seen in the panel. `NA` when too few samples cover the pair |
| `range_min_n_motif_kmers`, `range_max_n_motif_kmers` | median number of motif k-mers\* across samples at that length. `NA` when unavailable |

<p class="table-note">* An approximate measure of how much of the allele is made of the repeat motif.</p>

## Inputs for building your own catalog

Needed only if you decide to use your own tandem repeat list with `trcall prep-catalog` instead of using the pre-built catalog above.

### GRCh38 uniqmers

**hs38.ukm579.s02.sorted.txt.gz**

Every 31-mer that is conserved and unique in GRCh38 (i.e., uniqmer), with its position. Passed to `trcall prep-catalog -k`

| Column | Description |
| --- | --- |
| `chr` | chromosome |
| `position` | 0-based reference coordinate of the 31-mer |
| `kmer_id` | uniqmer ID: `a`, `x`, or `y` followed by numbers, for autosomes, chrX, and chrY |
| `strand` | reference strand, `+` or `-` |

### Uniqmer ID to sequence map

**ukm579.s02.txt.gz**

The sequence behind each uniqmer ID. Passed to `trcall prep-catalog -s`.

| Column | Description |
| --- | --- |
| `kmer_id` | uniqmer ID, as in the uniqmers file |
| `kmer_seq` | the 31-mer sequence as the reference lists it |