#!/usr/bin/env nextflow

process PROCESS_CACHE {
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
  cache_dir = meta.cache_dir
  force_create_config = params.force_create_config ? "--force" : ""
  
  '''
  process_cache.py \
    !{species} \
    !{assembly} \
    !{version} \
    --ini_file !{ini_file} \
    --cache_dir !{cache_dir} \
    !{force_create_config}
  '''
}
