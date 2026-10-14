# Install

trcall can be installed through Bioconda or built from source.

## Conda install

<!-- TODO: remove this admonition once the bioconda recipe is merged -->
!!! warning "Not yet available"
    The Bioconda package is not published yet. Build from source in the meantime.

```sh
conda install -c conda-forge -c bioconda trcall
```

The conda package also installs `samtools`. Other tools that are required for certain inputs are listed under [External programs](#external-programs).

## Build from source

**Requirements**:

- A C++17 compiler: GCC 8+, Clang 7+, or Apple Clang (Xcode 11+)
- zlib, with development headers (`zlib1g-dev` on Debian/Ubuntu, `zlib-devel` on
  RHEL/Fedora, included with Xcode on macOS)
- GNU make
- Linux or macOS (Windows through WSL)

```sh
git clone https://github.com/erkim11/trcall.git
cd trcall
make
```

Then follow the [Quick start](quick-start.md) to check that it runs.

## External programs

These must be on your PATH, but only for the inputs listed. None is required to build `trcall` or to run it on local SAM, BAM, FASTA, or FASTQ files.

| Tool | Needed for |
| --- | --- |
| `samtools` | CRAM input, alignments streamed from a URL (s3://, https://, gs://, ...) |
| `curl` | FASTA/FASTQ fetched over http, https, or ftp |
| `aws` CLI | FASTA/FASTQ fetched from s3:// |
| `gcloud` | FASTA/FASTQ fetched from gs:// |
| `gzip` | `trcall allele` output whose file name ends in `.gz` |
