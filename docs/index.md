# Introduction

TRumut is an alignment-free tandem repeat caller for long reads.

- Please refer to our preprint for why we built TRumut, how it works, and what it finds.
- Source code is available on [GitHub](https://github.com/erkim11/trumut).

## Workflow

![TRumut workflow](workflow.svg#only-light)
![TRumut workflow](workflow_dark_docs.svg#only-dark)

**Main workflow**:

Run `trumut allele` on each sample. Then, take its output as the input to `trumut cancer`, `trio`, or `use-pon`.

- `trumut cancer`: reports alleles in a tumor that are absent from its matched normal
- `trumut trio`: reports alleles in a child that are absent from both parents
- `trumut use-pon`: reports alleles outside the range seen in a panel of normals. Use it when there is no matched sample to compare against: a tumor without a matched normal, or a single tissue sample screened for mosaic variants

**Not needed for the main workflow**:

Usage and output details for these three are under [Advanced info](advanced/prep-catalog.md).

- `trumut prep-catalog`: builds a uniqmer catalog for `trumut allele` from your own tandem repeat catalog  
  &emsp;→ ready-made version of the catalog available under [Pre-built files](pre-built-files.md)
- `trumut build-pon`: builds a panel of normals for `trumut use-pon` from your own samples  
  &emsp;→ ready-made version of the panel available under [Pre-built files](pre-built-files.md)
- `trumut mutual`: compares two or more samples symmetrically, reporting alleles present in one sample but absent from another. The manuscript uses it to profile HiFi error rates

## Citation

TO-DO: add citation here.
