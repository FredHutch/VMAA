#!/usr/bin/env nextflow

// Container versions
container__cutadapt = "quay.io/fhcrc-microbiome/cutadapt:cutadapt_2.3_bcw_0.3.1"
container__bwa = "quay.io/fhcrc-microbiome/bwa:bwa.0.7.17__bcw.0.3.0I"
container__assembler = "quay.io/biocontainers/megahit:1.2.9--h8b12597_0"
container__ubuntu = "ubuntu:20.04"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
container__kraken2 = "quay.io/fhcrc-microbiome/kraken2:v2.0.9-beta"

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
  container "${container__ubuntu}"
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
    tuple val(sample_name), file(FASTQ)

    output:
    tuple val(sample_name), file("${sample_name}.cutadapt.fq.gz") optional true

"""
set -e 

cutadapt \
--minimum-length ${params.min_len} \
-j ${task.cpus} \
--cut ${params.bases_before_adapter} \
-g ^${params.adapter_F}...${params.adapter_R} \
-q ${params.qual_threshold},${params.qual_threshold} \
--discard-untrimmed \
--max-n ${params.max_n_prop} \
--times 3 \
--trim-n \
${FASTQ} | \
awk '{if(NR % 4 != 0){printf \$1 "\t"}else{print \$1}}' | \
grep -v ${params.adapter_F} | \
grep -v ${params.adapter_R} | \
tr '\t' '\n' | gzip -c > \
${sample_name}.cutadapt.fq.gz

# If the file is empty, delete it entirely
if (( \$(gunzip -c ${sample_name}.cutadapt.fq.gz | wc -l) == 0 )); then
  rm ${sample_name}.cutadapt.fq.gz
fi
"""
}

// Process to remove human reads
process remove_human {
    tag "Remove human reads"
    container "${container__bwa}"
    errorStrategy 'retry'
    label 'mem_veryhigh'


    input:
        tuple val(sample_name), file("${sample_name}.input.fq.gz")
        file hg_index_tgz

    output:
        tuple val(sample_name), file("${sample_name}.fq.gz") optional true

"""
#!/bin/bash

set -e

# Get the index name from the first file in the tar
bwa_index_prefix=\$(tar -ztvf ${hg_index_tgz} | head -1 | sed \'s/.* //\')

# Remove whichever ending this file has
bwa_index_prefix=\${bwa_index_prefix%.amb}
bwa_index_prefix=\${bwa_index_prefix%.ann}
bwa_index_prefix=\${bwa_index_prefix%.bwt}
bwa_index_prefix=\${bwa_index_prefix%.fai}
bwa_index_prefix=\${bwa_index_prefix%.pac}
bwa_index_prefix=\${bwa_index_prefix%.sa}

echo BWA index prefix is \${bwa_index_prefix}
echo "If this index prefix is not correct, consider remaking the tarball"
echo "so that it doesn't include anything other than the index files"

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
${sample_name}.input.fq.gz

echo Checking if alignment is empty  
if [[ -s alignment.sam ]]; then

  echo Extracting Unaligned Pairs 
  samtools \
    fastq \
    alignment.sam \
    --threads ${task.cpus} \
    -f 4 \
    -n \
    | gzip -c > ${sample_name}.fq.gz

  
  echo "Checking to see how many reads pass the human filtering"

  gunzip -c ${sample_name}.fq.gz | wc -l

  if (( \$(gunzip -c ${sample_name}.fq.gz | wc -l) == 0 )); then

    echo "Removing empty output file"
    rm ${sample_name}.fq.gz

  fi

else

  echo "Alignment SAM file was empty"

fi

echo Done 
"""
}

// De novo assembly
process assemble {
    tag "De novo metagenomic assembly"
    container "${container__assembler}"
    label 'mem_veryhigh'
    errorStrategy "retry"
    publishDir params.output_folder

    input:
        path fastq_list
    
    output:
        file "${params.output_prefix}.fasta.gz"
    
"""
set -e 

echo -e "Concatenating inputs"
cat ${fastq_list} > INPUT.fastq.gz

date
echo -e "Running Megahit\\n"

megahit \
    -r INPUT.fastq.gz \
    -o OUTPUT \
    -t ${task.cpus} \
    --k-min ${params.k_min} \
    --k-max ${params.k_max} \
    --k-step ${params.k_step} \
    --min-contig-len ${params.min_len} \

date
echo -e "\\nMaking sure output files are not empty\\n"
[[ \$(cat OUTPUT/final.contigs.fa | wc -l) > 0 ]]

date
echo -e "\\nRenaming output files\\n"

# Rename the output file
cat OUTPUT/final.contigs.fa | gzip -c > ${params.output_prefix}.fasta.gz

date
echo -e "\\nDone\\n"
"""
}

// Count the number of reads
process countReads {
  container "${container__ubuntu}"
  label 'io_limited'
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
  container "${container__ubuntu}"
  label 'io_limited'
  errorStrategy 'retry'
  publishDir params.output_folder
  
  input:
  file csv_list
  
  output:
  file "${params.output_prefix}.counts.csv"

  """
#!/bin/bash

set -e

echo "specimen,nreads,label" > ${params.output_prefix}.counts.csv
cat ${csv_list} >> ${params.output_prefix}.counts.csv
  """

}

process index {
  container "${container__bwa}"
  label "mem_medium"
  errorStrategy 'retry'

  input:
  file genome_fasta
  
  output:
  file "${genome_fasta}.tar"
  
  """
#!/bin/bash
set -e

echo "Indexing ${genome_fasta}"
bwa index ${genome_fasta}

echo "Tarring up indexed genome"
tar cvf ${genome_fasta}.tar ${genome_fasta}*
echo "Done"
  """

}

process align {
  container "${container__bwa}"
  label "mem_veryhigh"
  errorStrategy 'retry'
  publishDir "${params.output_folder}/bam/", mode: 'copy'

  input:
  tuple val(sample_name), file(input_fastq)
  file genome_index
  
  output:
  tuple val(sample_name), file("${sample_name}.bam")

  """
#!/bin/bash
set -e

echo Unpacking ${genome_index}
tar xvf ${genome_index}

echo "Running BWA"
bwa mem \
  -t ${task.cpus} \
  ${genome_index.name.replaceAll(/.tar/, "")} \
  ${input_fastq} \
  | samtools view -b -F 4 -o ${sample_name}.bam

# Sort the BAM file
samtools sort -o ${sample_name}.SORTED ${sample_name}.bam
mv ${sample_name}.SORTED.bam ${sample_name}.bam

# Index the BAM file
samtools index ${sample_name}.bam

echo Done aligning ${sample_name}.bam

    """

}

process calcStats {
  container "${container__bwa}"
  label "io_limited"
  errorStrategy 'retry'

  input:
  tuple val(sample_name), file(bam)
  
  output:
  tuple val(sample_name), file("${sample_name}.idxstats"), file("${sample_name}.stats"), file("${sample_name}.pileup"), file("${sample_name}.positions")

  """
#!/bin/bash
set -e

samtools sort ${bam} > ${sample_name}.sorted
samtools stats ${sample_name}.sorted > ${sample_name}.stats
samtools index ${sample_name}.sorted
samtools idxstats ${sample_name}.sorted > ${sample_name}.idxstats
samtools mpileup ${sample_name}.sorted > ${sample_name}.pileup

# Make a file with four columns, the contig, the bitwise flag, and the leftmost position, and the length of the mapped segment
samtools view ${sample_name}.sorted | awk '{print \$3 "\\t" \$2 "\\t" \$4 "\\t" length(\$10)}' > ${sample_name}.positions

rm ${sample_name}.sorted

  """

}

process summarize {
  container "${container__pandas}"
  label "mem_medium"
  errorStrategy 'retry'

  input:
  tuple val(sample_name), file(idxstats), file(stats), file(pileup), file(positions)
  
  output:
  file "${sample_name}.csv.gz"

  """
#!/usr/bin/env python3
import logging
import os
import gzip
import json
import pandas as pd
from math import log as ln

logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [nf-viral-metagenomics] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

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
logging.info("Reading in ${positions}")
positions = pd.read_csv(
  "${positions}", 
  sep="\\t", 
  header=None,
  names = [
    "contig", "flags", "pos", "len"
  ]
)

# Parse the flags columns from the positions file
logging.info("Parsing flags in ${positions}")
positions = pd.concat([positions, pd.DataFrame(map(parse_flags, positions["flags"]))], axis=1)

# Find the read start position
# If the read is aligned in the forward direction, use the leftmost position, otherwise use the rightmost
positions = positions.assign(
  start_pos = positions.apply(
    lambda r: r["pos"] + r["len"] if r["rc"] else r["pos"],
    axis=1
  )
)

# Calculate Shannon diversity per contig
logging.info("Calculating Shannon diversity per contig")
sdi = positions.groupby(
  "contig"
).apply(
  lambda contig_positions: shannon_divesity(contig_positions["start_pos"].value_counts().values)
)

# Read in the pileup
logging.info("Reading in ${pileup}")
pileup = pd.read_csv(
  "${pileup}", 
  sep="\\t", 
  header=None,
  names = ["contig", "position", "base", "depth", "aligned_bases", "aligned_quality"]
)

# Get the stats for each contig
logging.info("Reading in ${idxstats}")
contig_stats = pd.read_csv(
  "${idxstats}", 
  sep="\\t", 
  header=None,
  names = ["contig", "len", "mapped", "unmapped"]
).set_index(
  "contig"
)

# Number of positions with any aligned reads per contig
logging.info("Number of positions with any aligned reads per contig")
covlen = pileup.query(
  "depth > 0"
)["contig"].value_counts()

# Error rate over the entire alignment
logging.info("Reading overall error rate")
nerror = float(read_line("${stats}", "SN\\terror rate:").split("\\t")[0])

# Number of bases aligned per contig
nbases = pileup.groupby(
  "contig"
).apply(
  lambda df: df["depth"].sum()
)

# Make a single output object
logging.info("Making an output object")
output = dict([
  ("specimen", "${sample_name}"),
  ("reads_aligned", positions["contig"].value_counts()),
  ("bases_aligned", nbases),
  ("bases_covered", covlen),
  ("contig_length", contig_stats["len"]),
  ("error", nerror),
  ("entropy", sdi),
  ("specimen_total_reads", contig_stats.reindex(columns=["mapped", "unmapped"]).sum().sum())
])

# Format as a DataFrame
logging.info("Formatting as DataFrame")
output = pd.DataFrame(
  output
).drop(  # Drop the pseudo contig name used for unaligned reads
  index = "*"
)

# Add more summary columns, and reset the index used for the contig name
logging.info("Adding summary stats")
output = output.assign(
  depth = output["bases_aligned"] / output["contig_length"],
  coverage = output["bases_covered"] / output["contig_length"],
  proportion_of_reads = output["reads_aligned"] / output["specimen_total_reads"]
).reset_index(
).rename(
  columns = dict([("index", "contig")])
)

# Write out as CSV
logging.info("Writing out to ${sample_name}.csv.gz")
output.to_csv(
  "${sample_name}.csv.gz", 
  index=None
)

  """

}

process collect {
  container "${container__pandas}"
  label "io_limited"
  errorStrategy 'retry'
  publishDir params.output_folder

  input:
  file summary_csv_list

  output:
  file "${params.output_prefix}.summary.csv.gz"

  """
#!/usr/bin/env python3
import logging
import os
import gzip
import json
import pandas as pd

logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [nf-viral-metagenomics] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

summary_csv_list = "${summary_csv_list}".split(" ")
logging.info("Reading in a list of %d summary CSV files" % len(summary_csv_list))
summary_df = pd.concat([
  pd.read_csv(fp)
  for fp in summary_csv_list
])

logging.info("Saving to ${params.output_prefix}.summary.csv.gz")
summary_df.to_csv("${params.output_prefix}.summary.csv.gz", index=None)

logging.info("Done")

"""

}


process collect_with_kraken {
  container "${container__pandas}"
  label "io_limited"
  errorStrategy 'retry'
  publishDir params.output_folder

  input:
  file summary_csv_list
  file kraken2_tsv

  output:
  file "${params.output_prefix}.summary.csv.gz"

  """
#!/usr/bin/env python3
import logging
import os
import gzip
import json
import pandas as pd

logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [nf-viral-metagenomics] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

summary_csv_list = "${summary_csv_list}".split(" ")
logging.info("Reading in a list of %d summary CSV files" % len(summary_csv_list))
summary_df = pd.concat([
  pd.read_csv(fp)
  for fp in summary_csv_list
])

# Read in the Kraken2 results
kraken_df = pd.read_csv(
  "${kraken2_tsv}", 
  sep="\\t", 
  header=None,
  names = [
    "success",
    "query",
    "tax_id",
    "query_len",
    "hit_string"
  ]
).query(
  "success == 'C'"
).drop(
  columns=["success", "query_len"]
).set_index(
  "query"
)

# Add the Kraken results to the summary DataFrame
summary_df = summary_df.assign(
  tax_id = summary_df["contig"].apply(
    kraken_df["tax_id"].get
  ).fillna(
    0
  ).apply(
    int
  )
)
summary_df = summary_df.assign(
  hit_string = summary_df["contig"].apply(
    kraken_df["hit_string"].get
  )
)

logging.info("Saving to ${params.output_prefix}.summary.csv.gz")
summary_df.to_csv("${params.output_prefix}.summary.csv.gz", index=None)

logging.info("Done")

"""

}


// Perform taxonomic classification with Kraken2
process kraken2 {
  container "${container__kraken2}"
  errorStrategy 'retry'
  label "mem_veryhigh"

  input:
  file input_fasta
  file "DB/*"
  
  output:
  file "${params.output_prefix}.kraken2.gz"

  """
#!/bin/bash

set -e

ls -lahtr

ls -lahtr DB/

echo "Running kraken2"

KRAKEN2_DB_PATH=\$PWD kraken2 \
    --db DB \
    --threads ${task.cpus} \
    ${input_fasta} \
    | gzip -c > ${params.output_prefix}.kraken2.gz

  """

}