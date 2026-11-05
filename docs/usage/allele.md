# trumut allele

A worked example on simulated data is under [Quick start](../quick-start.md).

## Basic command

```sh
trumut allele [reads.bam] [catalog.tsv] [output.tsv] --input-type reads
trumut allele [assembly.fasta] [catalog.tsv] [output.tsv] --input-type assembly --ploidy [2]
```

- `input (reads.bam/assembly.fasta)`:

    - accepted types:

        - alignments (SAM, SAM.gz, BAM, CRAM)  
        - assemblies (FASTA)  
        - unaligned reads (FASTQ, FASTA, unaligned SAM/BAM)  

    - can be given as local file, stdin, or a URL

- `[catalog.tsv]`: uniqmer-annotated tandem repeat catalog from [Pre-built files](../pre-built-files.md) or `trumut prep-catalog`
- `[output.tsv]`: output path
- <code>--input-type reads&#124;assembly</code>: `reads` for sequencing reads. `assembly` for assembled contigs
- `--ploidy N`: required with `assembly`. Number of haplotype contigs per locus (e.g., `1` for a haploid assembly, `2` for diploid)

## More options

### For input

| Argument | Default | Description |
| --- | --- | --- |
| `--input-fmt STR` | sniffed | specify input file format |
| `-r, --reference FILE` |  | reference FASTA, needed for CRAM |

### For output

| Argument | Default | Description |
| --- | --- | --- |
| `--save-allele-seq` | off | add an `allele_seqs` column: the read sequence between the two uniqmers |
| `--save-query-names` | off | add a `query_names` column: read names, or contig names for assembly input |
| `--save-kmer-coords` | off | add `kmerA_ref_coords` and `kmerB_ref_coords` columns: the reference-genome position of each uniqmer |
| `-w, --wrong-kmer-hits FILE` | none | uniqmer pairs found in a read / contig on unexpected strands or in the wrong order. |
| `-m, --many-kmer-hits FILE` | none | uniqmer pairs that were found more than once in a single read / contig |
| `-d, --debug-file FILE` | none | one row per catalog uniqmer pair: whether it was dropped from the main output and, if so, why |

### Per input type

#### Indexed alignments

→ BAM, CRAM, or URL with a remote index

| Argument | Default | Description |
| --- | --- | --- |
| `-t, --target STR` | all | process one chromosome only, spelled exactly as in the BAM/CRAM header |
| `-j, --jobs INT` | 3 | number of chromosomes processed in parallel, one worker each; lowered automatically if that many would not fit in memory |
| `-q, --min-mapq INT` | 0 | skip alignments below this MAPQ |
| `-@, --threads INT` | 3 | threads passed to samtools for CRAM decoding, per worker |

#### Unindexed streams

→ SAM (aligned or not), unindexed BAM, stdin, FASTQ/FASTA reads

| Argument | Default | Description |
| --- | --- | --- |
| `-t, --target STR` | all | SAM/BAM streams only: process one chromosome. `chr1` and `1` both match |
| `-j, --jobs INT` | 3 | scanning threads |
| `-q, --min-mapq INT` | 0 | SAM/BAM streams only: skip alignments below this MAPQ |
| `-@, --threads INT` | 3 | BGZF decompression threads for bgzip FASTQ or unaligned BAM |

#### Assemblies

→ FASTA contigs

| Argument | Default | Description |
| --- | --- | --- |
| `-t, --target STR` | all | process one record (contig) only, by name |
| `-@, --threads INT` | 3 | BGZF decompression threads for a bgzip assembly |
