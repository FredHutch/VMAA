# Nextflow tools for viral metagenomics


### Goal

To create a reproducible pipeline for viral metagenomics that is
relatively portable between computational infrastructures.


### Foundation

Nextflow, a workflow language that enables researchers to run arbitrary
BASH commands within Docker images, creating reproducible workflows
which can be executed on any computational system, supporting local 
testing and deployment on an HPC or in the cloud.


### Data Storage

We will assume that all of the input data is being stored in the
Synapse object storage system. The workflow may be expanded in the 
future to support a larger collection of storage options.


### Usage

The list of samples for a given experiment will be stored in a text
file, one Synapse UUID per line.

The `viral-metagenomics.nf` workflow will then be executed on that
set of samples with the following command.

```
nextflow run --sample_sheet sample_sheet.txt viral-metagenomics.nf
```

By default execution will take place locally, but adding a configuration
file for Nextflow will then enable execution on cloud or HPC resources.
