# Nextflow tools for viral metagenomics


### Goal

To create a reproducible pipeline for viral metagenomics that is
relatively portable between computational infrastructures.


### Foundation

Nextflow, a workflow language that enables researchers to run arbitrary
BASH commands within Docker images, creating reproducible workflows
which can be executed on any computational system, supporting local 
testing and deployment on an HPC or in the cloud.


### Analysis

Quite simply, all input files will be aligned against all reference
genome sequences using the BWA algorithm (which is appropriate for
short-read metagenomics). A summary of the alignment against each 
genome will be written to a single CSV, containing information on the 
depth and coverage of sequencing.


### Data Storage

We will assume that all of the input data is being stored in a single
folder, either in local storage or in economy cloud storage (AWS S3).
In either case, all of the FASTQ files from a particular folder will
be analyzed against a set of reference genomes.

### Reference Genomes

The only requirement for the reference genomes is that you provide a
CSV with a column labeled 'accession' which corresponds to a NCBI nucleotide
accession. That reference will be downloaded and used for alignment.


### Output File Format

A summary of the alignment of all samples against all references will be
written to a file, whose name can be specified with the `--output_csv` flag.
That file will be placed in the `--output_directory` folder.


### Usage

```
nextflow \
    run \
        viral-metagenomics.nf \
        --input_directory <INPUT_DIRECTORY> \
        --output_directory <OUTPUT_DIRECTORY> \
        --viral_genome_csv <VIRAL_GENOME_CSV> \
        --output_csv <NAME_OF_OUTPUT_CSV>
```
