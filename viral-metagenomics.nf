#!/usr/bin/env nextflow

sample_sheet_ch = Channel.from(file(params.sample_sheet).readLines())
synapse_username = params.synapse_username
synapse_password = params.synapse_password

process fetch_from_synapse {

  container "quay.io/biocontainers/synapseclient@sha256:fc96a0c4cf72ff143419314e8b23cea1d266d495c55f45b1901fa0cc77e67153"
  cpus 1
  memory "2 GB"

  input:
  val synapse_uuid from sample_sheet_ch

  output:
  file "${synapse_uuid}.fastq.gz" into sample_fastq_ch

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

process make_kraken_db {

  cpus 32
  memory "240 GB"
  container "quay.io/fhcrc-microbiome/kraken2@sha256:e71b312bd428174b93377a8aaaec72a77524aba8b138cf569b929c44affee0c9"
  cpus 16
  memory "120 GB"

  output:
  file "KRAKEN_DB.tar"

  """
  kraken2-build --standard --threads 16 --db KRAKEN_DB && \
  tar cvf KRAKEN_DB.tar KRAKENDB/
  """

}
