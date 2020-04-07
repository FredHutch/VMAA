#!/usr/bin/env nextflow

// Preprocessing options
params.adapter_F = "CTGTCTCTTATACACATCT"
params.adapter_R = "CTGTCTCTTATACACATCT"

// Container versions
container__cutadapt = "quay.io/fhcrc-microbiome/cutadapt:cutadapt_2.3_bcw_0.3.1"
container__bwa = "quay.io/fhcrc-microbiome/bwa:bwa.0.7.17__bcw.0.3.0I"
container__assembler = "quay.io/biocontainers/megahit:1.2.9--h8b12597_0"

// Combine files which share a specimen label
workflow join_fastqs_by_specimen {
    take:
        fastq_ch

    main:
        fastq_ch.map {
            r -> [r["specimen"], file(r["fastq"])]
        }.groupTuple(
        ).branch {  // Split up the samples which have multiple FASTQ files
            single: it[1].size() == 1
            multiple: it[1].size() > 1
        }.set {
            grouped_fastq
        }

        // Join the FASTQ files for those samples
        joinFastq(
            grouped_fastq.multiple
        )

    emit:
        // Return a channel with every specimen represented as a tuple
        // with the specimen name and the FASTQ file
        joinFastq.out.mix(
          grouped_fastq.single.map {
            it -> [it[0], it[1][0]]
          }
        )
}

// Count the number of reads
process joinFastq {
  container "ubuntu:16.04"
  errorStrategy 'retry'
  
  input:
  tuple val(specimen), path(fastq_list)
  
  output:
  tuple val(specimen), path("${specimen}.fastq.gz")

  """
#!/bin/bash

set -e

cat ${fastq_list} > TEMP && mv TEMP "${specimen}.fastq.gz"
  """

}

// Process to run catadapt
process cutadapt {
    tag "Trim adapters from WGS reads"
    container "${container__cutadapt}"
    label 'mem_medium'
    errorStrategy 'retry'
    maxRetries 10

    input:
    tuple sample_name, file(FASTQ)

    output:
    tuple sample_name, file("${sample_name}.cutadapt.fq.gz")

"""
set -e 

cutadapt \
-j ${task.cpus} \
-a ${params.adapter} \
-q ${params.qual_threshold},${params.qual_threshold} \
-m ${params.min_len} \
--max-n ${params.max_n_prop} \
--trim-n \
-o ${sample_name}.cutadapt.fq.gz \
${FASTQ}
"""
}


// Process to remove human reads
process remove_human {
    tag "Remove human reads"
    container "${container__bwa}"
    errorStrategy 'retry'
    label 'mem_veryhigh'


    input:
        tuple sample_name, file(FASTQ)
        file hg_index_tgz

    output:
        tuple sample_name, file("${sample_name}.nohuman.fq.gz")

"""
#!/bin/bash

set -e

bwa_index_fn=\$(tar -ztvf ${hg_index_tgz} | head -1 | sed \'s/.* //\')
bwa_index_prefix=\${bwa_index_fn%.*}

echo BWA index prefix is \${bwa_index_prefix}

echo Extracting BWA index

mkdir -p hg_index/ 

tar -I pigz -xf ${hg_index_tgz} -C hg_index/

echo Files in index directory: 
ls -l -h hg_index 

# Make sure that there are files in the index directory
(( \$(ls -l -h hg_index | wc -l) > 0 ))

echo Running BWA 

bwa mem -t ${task.cpus} \
-T ${params.min_hg_align_score} \
-o alignment.sam \
hg_index/\$bwa_index_prefix \
${FASTQ}

echo Checking if alignment is empty  
[[ -s alignment.sam ]]
echo Extracting Unaligned Pairs 
samtools \
  fastq \
  alignment.sam \
  --threads ${task.cpus} \
  -f 4 \
  -n \
  | gzip -c > ${sample_name}.nohuman.fq.gz

echo Checking to see how many reads pass the human filtering
gunzip -c ${sample_name}.nohuman.fq.gz | wc -l
(( \$(gunzip -c ${sample_name}.nohuman.fq.gz | wc -l) > 0 ))

echo Done 
"""
}

// De novo assembly
process assemble {
    tag "De novo metagenomic assembly"
    container "${container__assembler}"
    label 'mem_veryhigh'
    errorStrategy "retry"

    input:
        path fastq_list
    
    output:
        file "contigs.fasta.gz"
    
"""
set -e 

echo -e "Concatenating inputs"
cat ${fastq_list} > INPUT.fastq.gz

date
echo -e "Running Megahit\\n"

megahit \
    -r INPUT.fastq.gz \
    -o OUTPUT \
    -t ${task.cpus}

date
echo -e "\\nMaking sure output files are not empty\\n"
[[ \$(cat OUTPUT/final.contigs.fa | wc -l) > 0 ]]

date
echo -e "\\nRenaming output files\\n"

# Rename the output file
cat OUTPUT/final.contigs.fa | gzip -c > contigs.fasta.gz

date
echo -e "\\nDone\\n"
"""
}

// Count the number of reads
process countReads {
  container "ubuntu:16.04"
  errorStrategy 'retry'
  
  input:
  tuple val(sample_name), file(fastq)
  
  output:
  file "${sample_name}.${params.count_reads_label}.counts.csv"

  """
#!/bin/bash

set -e

echo "Counting reads"
n=\$(gunzip -c "${fastq}" | awk 'NR % 4 == 1' | wc -l)

echo "Found \$n reads"
(( \$n > 0 ))

echo "${sample_name},\$n,${params.count_reads_label}" > "${sample_name}.${params.count_reads_label}.counts.csv"

  """

}

// Join the countReads CSV
process collectCountReads {
  container "ubuntu:16.04"
  errorStrategy 'retry'
  
  input:
  file csv_list
  
  output:
  file "counts.csv"

  """
#!/bin/bash

set -e

echo "specimen,nreads,label" > counts.csv
cat ${csv_list} >> counts.csv
  """

}

process index_viral_genomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  errorStrategy 'retry'

  input:
  file genome_fasta
  
  output:
  file "${genome_fasta}.tar"
  
  """
  bwa index ${genome_fasta}
  tar cvf ${genome_fasta}.tar ${genome_fasta}*
  """

}

process align_viral_genomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  errorStrategy 'retry'

  input:
  file genome_index
  each file(input_fastq)
  
  output:
  file "*.bam" optional true

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

  input:
  file bam
  
  output:
  set file("${bam}.idxstats"), file("${bam}.stats"), file("${bam}.pileup"), file("${bam}.positions")

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
  set file(idxstats), file(stats), file(pileup), file(positions)
  
  output:
  file "*json"

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
output["n_reads"] = n_reads
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
  file readcounts
  
  output:
  file "readcounts.csv"

  """
  echo file,n_reads > TEMP
  cat *csv >> TEMP && rm *csv && mv TEMP readcounts.csv
  """

}

process collect {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  errorStrategy 'retry'

  input:
  file all_jsons
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
