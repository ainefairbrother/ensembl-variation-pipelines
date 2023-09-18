#!/usr/bin/env nextflow

process PROCESS_CONSERVATION_DATA {
  cache false
  
  input:
  val meta

  output:
  val genome
  
  shell:
  genome = meta.genome
  version = params.version
  ini_file = params.ini_file
  conservation_data_dir = params.conservation_data_dir
  
  '''
  process_conservation_data.py \
    !{genome} \
    !{version} \
    --ini_file !{ini_file} \
    --conservation_data_dir !{conservation_data_dir}
  '''
}
