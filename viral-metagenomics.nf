#!/usr/bin/env nextflow

Channel.fromPath("${params.input_directory}*fastq.gz")
       .take(10)
       .into{ sample_fastq_kraken; sample_fastq_bwa; sample_fastq_viral_kraken}
kraken_db = file(params.kraken_db)
kraken_viral_db = file(params.kraken_viral_db)
viral_genome_ch = Channel.from(file(params.viral_genome_csv).readLines())
                         .splitCsv()
                         .map{it -> it[2]}
                         .filter( ~/^ftp.*/ )
                         .filter( ~/.*fna$/ )
                         .take(10)
                         .map{ it -> file(it) }

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


process viral_kraken {

  container "quay.io/fhcrc-microbiome/kraken2@sha256:ae4e647c2dd61c2f5595fd6682d50a4bde55fe9ba5e2ace424b70e68e205b8a6"
  cpus 32
  memory "240 GB"
  publishDir params.outdir
  errorStrategy 'retry'

  input:
  file kraken_viral_db
  file sample_fastq from sample_fastq_viral_kraken

  output:
  file "${sample_fastq}.viral.kraken.report.tsv"

  """
  set -e
  tar xvf ${kraken_viral_db}
  rm ${kraken_viral_db}

  kraken2 --db ${kraken_viral_db.simpleName} --threads 32 --report ${sample_fastq}.viral.kraken.report.tsv --output ${sample_fastq}.viral.kraken.tsv ${sample_fastq}

  rm -rf ${kraken_viral_db.simpleName}
  """
}

process index_viral_genomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 1
  memory "1 GB"
  errorStrategy 'retry'

  input:
  file genome_fasta from viral_genome_ch
  
  output:
  file "${genome_fasta}.tar" into indexed_genomes
  
  """
  bwa index ${genome_fasta}
  tar cvf ${genome_fasta}.tar ${genome_fasta}*
  """

}

process align_viral_genomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 8
  memory "8 GB"
  errorStrategy 'retry'

  input:
  file genome_index from indexed_genomes
  each file(input_fastq) from sample_fastq_bwa
  
  output:
  val record

  """
  tar xvf ${genome_index}
  genome_name=\$(echo ${genome_index} | sed 's/.tar//')
  bwa mem -t 8 \$genome_name ${input_fastq} | samtools view -b -F 4 -o ${input_fastq}.\$genome_name.bam
  """

}