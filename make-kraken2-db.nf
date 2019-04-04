#!/usr/bin/env nextflow

params.database_type = 'standard'

process make_kraken_db {

  // container "quay.io/fhcrc-microbiome/kraken2@sha256:ae4e647c2dd61c2f5595fd6682d50a4bde55fe9ba5e2ace424b70e68e205b8a6"
  container "job-definition://kraken2:1"
  cpus 32
  memory "220 GB"
  publishDir params.outdir
  scratch "/scratch"

  input:
  val dbname from params.dbname
  val database_type from params.database_type

  output:
  file "${dbname}.tar"

  """
  set -e;
  
  if [ "${database_type}" == "standard" ]
  then
    kraken2-build --standard --threads 32 --db ${dbname} --threads 32
  else
    echo "Downloading the taxonomy for ${dbname}"
    kraken2-build --download-taxonomy --db ${dbname} --threads 32
    echo ${database_type} | tr ',' '\\n' | while read domain
    do
      echo "Downloading library for \$domain"
      echo "kraken2-build --download-library \$domain --db ${dbname} --threads 32"
      kraken2-build --download-library \$domain --db ${dbname} --threads 32
    done
    echo "Building database"
    echo "kraken2-build --build --db ${dbname} --threads 32"
    kraken2-build --build --db ${dbname} --threads 32
  fi

  tar cvf ${dbname}.tar ${dbname}/*;
  rm -r ${dbname}/
  """
}
