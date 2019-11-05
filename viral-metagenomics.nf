#!/usr/bin/env nextflow

Channel.fromPath("${params.input_directory}*fastq.gz")
       .into{ count_reads; sample_fastq_bwa}

viral_accession_ch = Channel.from(file(params.viral_genome_csv))
                         .splitCsv(header: true)
                         .map{it -> it.accession}

process get_ncbi_accession{
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  errorStrategy 'retry'

  input:
  val accession from viral_accession_ch
  
  output:
  file "${accession}.fasta" into viral_genome_ch
  
  """
  wget -O ${accession}.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text"
  """
}

// Count the number of reads
process countReads {
  container "ubuntu:16.04"
  errorStrategy 'retry'
  
  input:
  file fastq from count_reads
  
  output:
  file "${fastq}.counts.csv" into counts

  afterScript "rm *"

  """
#!/bin/bash

set -e

gzip -t "${fastq}"

n=\$(gunzip -c "${fastq}" | awk 'NR % 4 == 1' | wc -l)

(( \$n > 0 ))

echo "${fastq},\$n" > "${fastq}.counts.csv"

  """

}

process index_viral_genomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
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
  errorStrategy 'retry'
  publishDir "${params.output_directory}/bam/"

  input:
  file genome_index from indexed_genomes.toSortedList().flatten().collate(10)
  each file(input_fastq) from sample_fastq_bwa
  
  output:
  file "*.bam" optional true into bam_ch

  afterScript "rm *"

  """
  for genome_index in *tar; do
    echo Processing \$genome_index
    tar xvf \$genome_index
    genome_name=\$(echo \$genome_index | sed 's/.tar//')
    bwa mem -t 8 \$genome_name ${input_fastq} | samtools view -b -F 4 -o ${input_fastq}.\$genome_name.bam
    echo Done aligning to \$genome_name
    
    # If zero reads were aligned, delete the BAM file
    [[ -f ${input_fastq}.\$genome_name.bam ]] && \
    [[ ! -s ${input_fastq}.\$genome_name.bam ]] && \
    rm ${input_fastq}.\$genome_name.bam
    echo Checked for empty file
    
    [[ -s ${input_fastq}.\$genome_name.bam ]] && \
    (( \$(samtools view ${input_fastq}.\$genome_name.bam | wc -l) == 0 )) && \
    rm ${input_fastq}.\$genome_name.bam
    echo Done with \$genome_name
  done

    """

}

process alignment_stats {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  errorStrategy 'retry'
  publishDir "${params.output_directory}/stats/"

  input:
  file bam from bam_ch.flatten()
  
  output:
  set file("${bam}.idxstats"), file("${bam}.stats"), file("${bam}.pileup"), file("${bam}.positions") into stats_ch

  """
  samtools sort ${bam} > ${bam}.sorted
  samtools stats ${bam}.sorted > ${bam}.stats
  samtools index ${bam}.sorted
  samtools idxstats ${bam}.sorted > ${bam}.idxstats
  samtools mpileup ${bam}.sorted > ${bam}.pileup

  # Make a file with three columns, the bitwise flag, and the leftmost position, and the length of the mapped segment
  samtools view ${bam}.sorted | awk '{print \$2 "\\t" \$4 "\\t" length(\$10)}' > ${bam}.positions

  rm ${bam}.sorted ${bam}

  """

}

process summarize_each {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  errorStrategy 'retry'

  input:
  set file(idxstats), file(stats), file(pileup), file(positions) from stats_ch
  
  output:
  file "*json" into all_stats_ch

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd
from math import log as ln

def read_line(fp, prefix):
    with open(fp, 'rt') as f:
        for line in f:
            if line.startswith(prefix):
                return line.replace(prefix, '').strip(" ").strip("\\t")

def parse_flags(int_flags):
    output = dict()
    for n, flag in [
        (2048, "supplementary"),
        (1024, "duplicate"),
        (512, "fail_filter"),
        (256, "secondary"),
        (128, "last"),
        (64, "first"),
        (32, "next_rc"),
        (16, "rc"),
        (8, "next_unmapped"),
        (4, "unmapped"),
        (2, "aligned_properly"),
        (1, "multiple_segments"),
    ]:
        if int_flags >= n:
            output[flag] = True
            int_flags = int_flags - n
        else:
            output[flag] = False
    assert int_flags == 0, int_flags
    return output


def shannon_divesity(counts):
    # See https://gist.github.com/audy/783125
    
    def p(n, N):
        # Relative abundance
        if n is 0:
            return 0
        else:
            return (float(n)/N) * ln(float(n)/N)
            
    N = sum(counts)
    
    return -sum(p(n, N) for n in counts if n is not 0)


# Read in the summary of alignment positions
positions = pd.read_csv("${positions}", sep="\\t", header=None)
positions = pd.concat([positions, pd.DataFrame(map(parse_flags, positions[0]))], axis=1)
# If the read is aligned in the forward direction, use the leftmost position, otherwise use the rightmost
position_list = positions.apply(
    lambda r: r[1] + r[2] if r["rc"] else r[1],
    axis=1
).tolist()

# Calculate Shannon diversity
sdi = shannon_divesity(pd.Series(position_list).value_counts().values)

pileup = pd.read_csv("${pileup}", sep="\\t", header=None)
n_reads = int(read_line("${stats}", "SN\\treads mapped:"))
reflen = int(pd.read_csv("${idxstats}", sep="\\t", header=None).iloc[0, 1])
covlen = pileup[3].shape[0]
nerror = float(read_line("${stats}", "SN\\terror rate:").split("\\t")[0])
nbases = pileup[3].sum()

output = dict()
output["file"] = "${pileup}".replace(".pileup", "")
output["depth"] = nbases / reflen
output["coverage"] = covlen / reflen
output["error"] = nerror
output["genome_length"] = reflen
output["n_reads"] = len(position_list)
output["entropy"] = sdi

json_fp = "${pileup}".replace(".pileup", ".json")
assert json_fp.endswith(".json")
with open(json_fp, "wt") as fo:
    fo.write(json.dumps(output, indent=4))

assert os.path.exists(json_fp)

  """

}

process collectCounts {
  container "ubuntu:16.04"
  errorStrategy 'retry'

  input:
  file readcounts from counts.collect()
  
  output:
  file "readcounts.csv" into readcounts_csv

  """
  echo file,n_reads > TEMP
  cat *csv >> TEMP && rm *csv && mv TEMP readcounts.csv
  """

}

process collect {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  publishDir "${params.output_directory}/"
  errorStrategy 'retry'

  input:
  file all_jsons from all_stats_ch.collect()
  file readcounts_csv
  
  output:
  file "${params.output_csv}"

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd

readcounts = pd.read_csv("${readcounts_csv}").set_index(
    "file"
)["n_reads"].to_dict()

def match_file_name(file_name):
    match = None
    for k, v in readcounts.items():
        if file_name.startswith(k):
            assert match is None, "Duplicate file name matching: %s" % (file_name)
            match = v
    return match

df = pd.DataFrame([
    json.load(open(fp, "rt"))
    for fp in os.listdir(".")
    if fp.endswith(".json")
])

df["total_reads"] = df["file"].apply(match_file_name)
assert df["total_reads"].isnull().sum() == 0
assert (df["total_reads"] > 0).all()
df["prop_reads"] = df["n_reads"] / df["total_reads"]

df.to_csv("${params.output_csv}", index=None)

"""

}
