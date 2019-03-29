#!/usr/bin/env nextflow

sample_sheet_ch = Channel.from(file(params.sample_sheet).readLines())
synapse_username = params.synapse_username
synapse_password = params.synapse_password
params.fve_db_idx_synapse = "syn18378932"
params.fve_db_list_synapse = "syn18378769"
kraken_db = file(params.kraken_db)
params.fve_cr = "0.1"
params.fve_co = "0.01"
params.fve_cn = "1"

process fetch_sample_fastq {

  container "quay.io/biocontainers/synapseclient@sha256:fc96a0c4cf72ff143419314e8b23cea1d266d495c55f45b1901fa0cc77e67153"
  cpus 1
  memory "2 GB"
  errorStrategy 'retry'

  input:
  val synapse_uuid from sample_sheet_ch

  output:
  file "${synapse_uuid}.fastq.gz" into sample_fastq_fve, sample_fastq_kraken

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
  errorStrategy 'retry'

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
  errorStrategy 'retry'

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

process fast_virome_explorer {

  container "quay.io/fhcrc-microbiome/fastviromeexplorer@sha256:555103371bc4b21be7fba64732e431f5bfc5ba2cf9305397ea8b4a5bb9a45f32"
  cpus 4
  memory "16 GB"
  errorStrategy 'retry'
  publishDir params.outdir

  input:
  file db_idx from fve_db_idx
  file db_list from fve_db_list
  file sample_fastq from sample_fastq_fve

  output:
  file "${sample_fastq}.fve.tsv"

  """
  set -e;
  
  java \
  -cp /usr/local/FastViromeExplorer/bin \
  FastViromeExplorer \
  -l ${db_list} \
  -1 ${sample_fastq} \
  -i ${db_idx} \
  -o ./ \
  -cr ${params.fve_cr} \
  -co ${params.fve_co} \
  -cn ${params.fve_cn}

  mv FastViromeExplorer-final-sorted-abundance.tsv ${sample_fastq}.fve.tsv;

  rm ${db_idx}

  """
}


process kraken {

  container "quay.io/fhcrc-microbiome/kraken2@sha256:ae4e647c2dd61c2f5595fd6682d50a4bde55fe9ba5e2ace424b70e68e205b8a6"
  cpus 32
  memory "240 GB"
  publishDir params.outdir
  errorStrategy 'retry'

  input:
  file kraken_db
  file sample_fastq from sample_fastq_kraken

  output:
  file "${sample_fastq}.kraken.report.tsv"

  """
  set -e
  tar xvf ${kraken_db}
  rm ${kraken_db}

  kraken2 --db ${kraken_db.simpleName} --threads 32 --report ${sample_fastq}.kraken.report.tsv --output ${sample_fastq}.kraken.tsv ${sample_fastq}

  rm -rf ${kraken_db.simpleName}
  """
}