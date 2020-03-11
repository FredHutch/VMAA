#!/usr/bin/env nextflow

// Using DSL-2
nextflow.preview.dsl=2

// Import modules used in the workflow
include join_fastqs_by_specimen from './modules/modules'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/VMAA <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)
      --output_folder       Folder to place analysis outputs
      --output_prefix       Text used as a prefix for output files

    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain the columns `specimen` and `fastq`. FASTQ files should be gzip-compressed.

    """.stripIndent()
}

workflow {

  // Show help message if the user specifies the --help flag at runtime
  params.help = false
  if (params.help || params.manifest == null || params.output_folder == null || params.output_prefix == null){
      // Invoke the function above which prints the help message
      helpMessage()
      // Exit out and do not run anything else
      exit 1
  }

  // Check to make sure that the manifest file can be found
  if (file(params.manifest).isEmpty()){
    log.info"""
    Cannot find ${params.manifest}
    """.stripIndent()
    exit 1
  }

  // Parse the manifest file
  fastq_ch = Channel.from(
    file(params.manifest).splitCsv(
          header: true, 
          sep: ","
      )
  )

  // First thing, join together files which are from the same specimen
  join_fastqs_by_specimen(
    fastq_ch
  )

  join_fastqs_by_specimen.out.view()

}