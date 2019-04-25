#!/usr/bin/env nextflow

sample_sheet_ch = Channel.from(file(params.sample_sheet).readLines())
synapse_username = params.synapse_username
synapse_password = params.synapse_password

process fetch_sample_fastq {

  container "quay.io/biocontainers/synapseclient@sha256:fc96a0c4cf72ff143419314e8b23cea1d266d495c55f45b1901fa0cc77e67153"
  cpus 1
  memory "2 GB"
  errorStrategy 'retry'
  publishDir "${params.outdir}"

  input:
  val synapse_uuid from sample_sheet_ch

  output:
  file "${synapse_uuid}.fastq.gz"

  """
  set -e;

  filename=\$(synapse -u ${synapse_username} -p ${synapse_password} get ${synapse_uuid} | grep Downloaded | sed 's/Downloaded file: //')
  echo "Downloaded \$filename"

  # Make sure the file exists
  [[ -s \$filename ]]

  gzip -t \$filename && \
  mv \$filename ${synapse_uuid}.fastq.gz || \
  gzip -c \$filename > ${synapse_uuid}.fastq.gz

  """  
}
