#!/usr/bin/env nextflow

sample_sheet_ch = Channel.from(file(params.sample_sheet).readLines())
synapse_username = params.synapse_username
synapse_password = params.synapse_password
params.fve_db_list_synapse = "syn18378769"

process fetch_sample_fastq {

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

process fetch_db_idx {

  container "quay.io/biocontainers/synapseclient@sha256:fc96a0c4cf72ff143419314e8b23cea1d266d495c55f45b1901fa0cc77e67153"
  cpus 1
  memory "2 GB"

  input:
  val synapse_uuid from params.fve_db_idx_synapse

  output:
  file "kallisto_db.idx" into fve_db_idx

  """
  set -e;

  filename=\$(synapse -u ${synapse_username} -p ${synapse_password} get ${synapse_uuid} | grep Downloaded | sed 's/Downloaded file: //')
  echo "Downloaded \$filename"
  mv \$filename kallisto_db.idx
  """  
}

process fetch_db_list {

  container "quay.io/biocontainers/synapseclient@sha256:fc96a0c4cf72ff143419314e8b23cea1d266d495c55f45b1901fa0cc77e67153"
  cpus 1
  memory "2 GB"

  input:
  val synapse_uuid from params.fve_db_list_synapse

  output:
  file "kallisto_db.list" into fve_db_list

  """
  set -e;

  filename=\$(synapse -u ${synapse_username} -p ${synapse_password} get ${synapse_uuid} | grep Downloaded | sed 's/Downloaded file: //')
  echo "Downloaded \$filename"
  mv \$filename kallisto_db.list
  """  
}

// process fast_virome_explorer {

//   // container "quay.io/fhcrc-microbiome/fastviromeexplorer@sha256:555103371bc4b21be7fba64732e431f5bfc5ba2cf9305397ea8b4a5bb9a45f32"
//   cpus 4
//   memory "16 GB"

//   input:
//   each file(sample_fastq) from sample_fastq_ch
//   file fve_db_ix
//   file fve_db_list

//   output:
//   "${sample_fastq}.fve.tsv"

//   """
//     java \
//     -cp /usr/local/FastViromeExplorer/bin \
//     FastViromeExplorer \
//     -l ${fve_db_list} \
//     -1 ${sample_fastq} \
//     -i ${fve_db} \
//     -o ./ && \
//     mv FastViromeExplorer-final-sorted-abundance.tsv ${sample_fastq}.fve.tsv
//   """
// }
