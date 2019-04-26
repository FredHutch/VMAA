#!/usr/bin/env nextflow

Channel.fromPath("${params.input_directory}*fastq.gz")
       .take(10)
       .into{ sample_fastq_kraken; sample_fastq_bwa; sample_fastq_viral_kraken}
kraken_db = file(params.kraken_db)
kraken_viral_db = file(params.kraken_viral_db)
viral_genome_ch = Channel.from(file(params.viral_genome_csv))
                         .splitCsv(header: true)
                         .map{it -> it.genome_ftp}
                         .take(10)
                         .map{ it -> file(it) }

// process kraken {

//   container "quay.io/fhcrc-microbiome/kraken2@sha256:ae4e647c2dd61c2f5595fd6682d50a4bde55fe9ba5e2ace424b70e68e205b8a6"
//   cpus 32
//   memory "240 GB"
//   publishDir "${params.output_directory}/kraken"
//   errorStrategy 'retry'

//   input:
//   file kraken_db
//   file sample_fastq from sample_fastq_kraken

//   output:
//   file "${sample_fastq}.kraken.report.tsv"

//   """
//   set -e
//   tar xvf ${kraken_db}
//   rm ${kraken_db}

//   kraken2 --db ${kraken_db.simpleName} --threads 32 --report ${sample_fastq}.kraken.report.tsv --output ${sample_fastq}.kraken.tsv ${sample_fastq}

//   rm -rf ${kraken_db.simpleName}
//   """
// }


// process viral_kraken {

//   container "quay.io/fhcrc-microbiome/kraken2@sha256:ae4e647c2dd61c2f5595fd6682d50a4bde55fe9ba5e2ace424b70e68e205b8a6"
//   cpus 32
//   memory "240 GB"
//   publishDir "${params.output_directory}/viral_kraken"
//   errorStrategy 'retry'

//   input:
//   file kraken_viral_db
//   file sample_fastq from sample_fastq_viral_kraken

//   output:
//   file "${sample_fastq}.viral.kraken.report.tsv"

//   """
//   set -e
//   tar xvf ${kraken_viral_db}
//   rm ${kraken_viral_db}

//   kraken2 --db ${kraken_viral_db.simpleName} --threads 32 --report ${sample_fastq}.viral.kraken.report.tsv --output ${sample_fastq}.viral.kraken.tsv ${sample_fastq}

//   rm -rf ${kraken_viral_db.simpleName}
//   """
// }

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
  publishDir "${params.output_directory}/bam/"

  input:
  file genome_index from indexed_genomes
  each file(input_fastq) from sample_fastq_bwa
  
  output:
  file "*.bam" optional true into bam_ch

  """
  tar xvf ${genome_index}
  genome_name=\$(echo ${genome_index} | sed 's/.tar//')
  bwa mem -t 8 \$genome_name ${input_fastq} | samtools view -b -F 4 -o ${input_fastq}.\$genome_name.bam
  echo Done aligning
  rm ${input_fastq}
  echo Removed input file
  
  # If zero reads were aligned, delete the BAM file
  [[ -f ${input_fastq}.\$genome_name.bam ]] && \
  [[ ! -s ${input_fastq}.\$genome_name.bam ]] && \
  rm ${input_fastq}.\$genome_name.bam
  echo Checked for empty file
  
  [[ -s ${input_fastq}.\$genome_name.bam ]] && \
  (( \$(samtools view ${input_fastq}.\$genome_name.bam | wc -l) == 0 )) && \
  rm ${input_fastq}.\$genome_name.bam
  echo Done
    """

}

process alignment_stats {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 1
  memory "1 GB"
  errorStrategy 'retry'
  publishDir "${params.output_directory}/stats/"

  input:
  file bam from bam_ch
  
  output:
  set file("${bam}.idxstats"), file("${bam}.stats"), file("${bam}.pileup") into stats_ch

  """
  samtools sort ${bam} > ${bam}.sorted
  samtools stats ${bam}.sorted > ${bam}.stats
  samtools index ${bam}.sorted
  samtools idxstats ${bam}.sorted > ${bam}.idxstats
  samtools mpileup ${bam}.sorted > ${bam}.pileup
  rm ${bam}.sorted ${bam}

  """

}

process summarize_each {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  cpus 1
  memory "4 GB"
  // errorStrategy 'retry'

  input:
  set file(idxstats), file(stats), file(pileup) from stats_ch
  
  output:
  file "*json" into all_stats_ch

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd

def read_line(fp, prefix):
    with open(fp, 'rt') as f:
        for line in f:
            if line.startswith(prefix):
                return line.replace(prefix, '').strip(" ").strip("\\t")

pileup = pd.read_csv("${pileup}", sep="\\t", header=None)
n_reads = int(read_line("${stats}", "SN\\treads mapped:"))
reflen = int(pd.read_csv("${idxstats}", sep="\\t", header=None).iloc[0, 1])
covlen = pileup[3].shape[0]
nerror = float(read_line("${stats}", "SN\\terror rate:").split("\\t")[0])
nbases = pileup[3].sum()

output = dict()
output["depth"] = nbases / reflen
output["coverage"] = covlen / reflen
output["error"] = nerror
output["genome_length"] = reflen
output["n_reads"] = n_reads

json_fp = "${pileup}".replace(".pileup", ".json")
assert json_fp.endswith(".json")
with open(json_fp, "wt") as fo:
    fo.write(json.dumps(output, indent=4))

assert os.path.exists(json_fp)

  """

}

process collect {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  cpus 1
  memory "4 GB"
  publishDir "${params.output_directory}/"
  // errorStrategy 'retry'

  input:
  file all_jsons from all_stats_ch.collect()
  
  output:
  file "${params.output_csv}"

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd

pd.DataFrame([
    json.load(open(fp, "rt"))
    for fp in os.listdir(".")
    if fp.endswith(".json")
]).to_csv("${params.output_csv}", index=None)
"""

}