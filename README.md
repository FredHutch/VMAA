# Viral Metagenomics by Nucleotide Alignment


### Goal

Reproducible computational pipeline for viral metagenomics that is
relatively portable between computational infrastructures.


### Approach

The analytical workflow will perform the following steps for each sample:

  1. Perform quality trimming and remove adapters
  2. Remove optical / PCR duplicates
  3. Remove human sequences by mapping against the human genome
  4. Assemble reads _de novo_ into contigs
  5. Map reads against contigs
  6. Calculate summary metrics for each contig
  7. Collect results across all samples
