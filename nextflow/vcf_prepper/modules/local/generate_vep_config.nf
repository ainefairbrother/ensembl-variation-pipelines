#!/usr/bin/env nextflow

process GENERATE_VEP_CONFIG {
  cache false
  
  input:
  val meta
    
  output:
  val genome
  
  shell:
  force_create_config = params.force_create_config
  genome = meta.genome
  version = params.version
  ini_file = params.ini_file
  vep_config = meta.vep_config
  cache_dir = params.cache_dir
  fasta = meta.fasta
  repo_dir = params.repo_dir
  
  '''
  if [[ ! -e !{vep_config} || !{force_create_config} == 1 ]]; then
    generate_vep_config.py \
      !{genome} \
      !{version} \
      --ini_file !{ini_file} \
      --vep_config !{vep_config} \
      --cache_dir !{cache_dir} \
      --fasta !{fasta} \
      --repo_dir !{repo_dir}
  fi
  '''
}
