#!/usr/bin/env nextflow

// Parameters used here
params.output_folder = false
params.output_prefix = false

// Container with kraken2 installed
container__kraken2 = "quay.io/fhcrc-microbiome/kraken2:v2.0.9-beta"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/build_kraken2_db <ARGUMENTS>
    
    Required Arguments:
      --output_folder       Folder to place database
      --output_prefix       Text used as a prefix for output files

    """.stripIndent()
}


// Build the Kraken2 database
process build_kraken2_db {
  container "${container__kraken2}"
  errorStrategy 'retry'
  publishDir "${params.output_folder}"
  memory 240.Gb
  cpus 32
  
  output:
  file "${params.output_prefix}*"

  """
#!/bin/bash

set -e

kraken2-build \
    --standard \
    --threads ${task.cpus} \
    --db ${params.output_prefix}

  """

}