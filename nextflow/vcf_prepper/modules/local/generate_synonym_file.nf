#!/usr/bin/env nextflow

process GENERATE_SYNONYM_FILE {
  cache false
  
  input:
  val meta
  
  output:
  val genome
  
  shell:
  genome = meta.genome
  species = meta.species
  assembly = meta.assembly
  version = params.version
  ini_file = params.ini_file
  synonym_file = meta.synonym_file
  force_create_config = params.force_create_config ? "--force" : ""
  
  '''
  generate_synonym_file.py \
    !{species} \
    !{assembly} \
    !{version} \
    --ini_file !{ini_file} \
    --synonym_file !{synonym_file} \
    !{force_create_config}
  '''
}
