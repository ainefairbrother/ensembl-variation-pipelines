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
  species = meta.species
  assembly = meta.assembly
  version = params.version
  ini_file = params.ini_file
  vep_config = meta.vep_config
  cache_dir = meta.cache_dir
  fasta_dir = meta.fasta_dir
  conservation_data_dir = meta.conservation_data_dir
  repo_dir = params.repo_dir
  
  '''
  if [[ ! -e !{vep_config} || !{force_create_config} == 1 ]]; then
    generate_vep_config.py \
      !{species} \
      !{assembly} \
      !{version} \
      --ini_file !{ini_file} \
      --vep_config !{vep_config} \
      --cache_dir !{cache_dir} \
      --fasta_dir !{fasta_dir} \
      --conservation_data_dir !{conservation_data_dir} \
      --repo_dir !{repo_dir}
  fi
  '''
}
