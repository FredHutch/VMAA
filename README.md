# Viral Metagenomics by Assembly and Alignment (VMAA)


### Goal

Reproducible computational pipeline for viral metagenomics that is
relatively portable between computational infrastructures.


### Approach

The analytical workflow will perform the following steps for each sample:

  1. Perform quality trimming and remove adapters
  2. Remove optical / PCR duplicates
  3. Remove human sequences by mapping against the human genome
  4. Assemble reads _de novo_ into contigs
  5. Map reads against contigs
  6. Calculate summary metrics for each contig
  7. Collect results across all samples


### Invoking the pipeline

The VMAA pipeline can be invoked by running `nextflow` against
this repository, e.g. `nextflow run FredHutch/VMAA <ARGS>`.

```
Usage:

nextflow run FredHutch/VMAA <ARGUMENTS>

Required Arguments:
  --manifest            CSV file listing samples (see below)
  --output_folder       Folder to place analysis outputs
  --output_prefix       Text used as a prefix for output files
  --human_genome_tar    Indexed human genome, gzipped tarball

Optional:
  --kraken2_folder      Folder containing Kraken2 database (omit trailing /)
  --kraken2_prefix      Prefix used to build Kraken2 database
  --k_min               Minimum k-mer size used for de novo assembly (default: 13)
  --k_max               Maximum k-mer size used for de novo assembly (default: 25)
  --k_step              Interval k-mer step size used for de novo assembly (default: 2)

Manifest:
  The manifest is a CSV with a header indicating which samples correspond to which files.
  The file must contain the columns `specimen` and `fastq`. FASTQ files should be gzip-compressed.
```
