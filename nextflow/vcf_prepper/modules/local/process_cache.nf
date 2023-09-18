#!/usr/bin/env nextflow

process PROCESS_CACHE {
  cache false
  
  input:
  val meta

  output:
  val genome
  
  shell:
  genome = meta.genome
  version = params.version
  ini_file = params.ini_file
  cache_dir = params.cache_dir
  
  '''
  process_cache.py \
    !{genome} \
    !{version} \
    --ini_file !{ini_file} \
    --cache_dir !{cache_dir}
  '''
}
