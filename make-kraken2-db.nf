#!/usr/bin/env nextflow

process make_kraken_db {

  // container "quay.io/fhcrc-microbiome/kraken2@sha256:ae4e647c2dd61c2f5595fd6682d50a4bde55fe9ba5e2ace424b70e68e205b8a6"
  container "job-definition://kraken2:1"
  cpus 32
  memory "220 GB"
  publishDir params.outdir
  scratch "/scratch"

  input:
  val dbname from params.dbname

  output:
  file "${dbname}.tar"

  """
  set -e;
  df -h;
  kraken2-build --standard --threads 32 --db ${dbname};
  tar cvf ${dbname}.tar ${dbname}/*;
  rm -r ${dbname}/
  """
}
