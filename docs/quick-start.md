# Quick start

Here we have a test run with small simulated samples so you can try `trumut` before using your own data. Follow the steps below to confirm everything runs.

## Get the data

First, download the test data. It lives in the test folder of the TRumut GitHub repository: [https://github.com/erkim11/trumut](https://github.com/erkim11/trumut)

```sh
# run this
git clone https://github.com/erkim11/trumut.git
cd trumut/test
```

## Run

!!! warning
    If you built TRumut from source, use `../trumut` in place of `trumut` in commands below.

```sh
# call tandem repeat alleles in each sample
trumut allele test_reads1.bam test_tr_catalog.tsv.gz reads1.alleles.tsv --input-type reads
trumut allele test_reads2.bam test_tr_catalog.tsv.gz reads2.alleles.tsv --input-type reads
trumut allele test_reads3.bam test_tr_catalog.tsv.gz reads3.alleles.tsv --input-type reads
```

```sh
# cancer: sample 3 plays the tumor, sample 1 its matched normal
trumut cancer reads3.alleles.tsv reads1.alleles.tsv output_cancer.tsv --input-type reads
```

```sh
# trio: sample 3 plays the child, samples 1 and 2 the parents
trumut trio reads3.alleles.tsv reads1.alleles.tsv reads2.alleles.tsv output_trio.tsv --input-type reads
```

```sh
# use-pon: sample 1 against a small panel of normals
trumut use-pon reads1.alleles.tsv test_pon.tsv output_pon.tsv --input-type reads
```

Each command should finish in seconds.

## Check

The first lines of each comparison output should look like this:

```sh
head -n 3 output_cancer.tsv | column -t
```

```
chr   region_start  region_end  kmerA_id     kmerB_id     loci                                                      allele_length  offset  n_motif_kmers_ref_median  n_motif_kmers_test_median  n_nonmotif_kmers_in  n_reads  n_reads_in_thres  event_type  phase_counts
chr1  776905        777048      a4799224467  a8869809255  776905-776921:AG,776913-777048:CCCCG                      168            11      112                       118                        nan                  5        12                clonal      0:7:5
chr2  986437        987067      a2892033592  a355770067   986437-987067:TCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTAT  734            -83     741.5                     663.5                      nan                  10       16                clonal      11:0:5
```

```sh
head -n 3 output_trio.tsv | column -t
```

```
chr   region_start  region_end  kmerA_id     kmerB_id    loci                                                      allele_length  offset  n_motif_kmers_ref_median  n_motif_kmers_test_median  n_nonmotif_kmers_in  n_reads  n_reads_in_thres  event_type  phase_counts
chr2  986437        987067      a2892033592  a355770067  986437-987067:TCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTAT  734            -82     739                       663.5                      nan                  10       16                clonal      11:0:5
chr2  986437        987067      a2892033592  a355770067  986437-987067:TCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTAT  608            40      500                       541                        nan                  6        11                clonal      0:11:0
```

```sh
head -n 3 output_pon.tsv | column -t
```

```
chr   region_start  region_end  kmerA_id     kmerB_id     loci                                  allele_length  offset  n_motif_kmers_ref_median  n_motif_kmers_test_median  n_reads  n_reads_in_thres  event_type      phase_counts
chr1  776905        777048      a4799224467  a8869809255  776905-776921:AG,776913-777048:CCCCG  187            13      119                       142                        6        12                clonal          8:0:4
chr2  1183385       1183470     a6387550638  a2329621088  1183385-1183470:CCCCG                 106            -10     73.5                      65                         3        20                private_phased  15:0:5
```

Each `reads*.alleles.tsv` file has 29 lines. If your output matches, the install is good to go. Column descriptions are under [Output](output/allele.md).

## Other input types

The same catalog works with unaligned reads and assemblies:

```sh
trumut allele test_reads1.fa.gz test_tr_catalog.tsv.gz reads1.fa.alleles.tsv --input-type reads
trumut allele test_assembly1.fa test_tr_catalog.tsv.gz asm1.alleles.tsv --input-type assembly --ploidy 2
```

Both give 29 lines.

## Next steps

For your own data, see [Usage](usage/allele.md) for each command and [Pre-built files](pre-built-files.md) for the genome-wide catalog and panel of normals.
