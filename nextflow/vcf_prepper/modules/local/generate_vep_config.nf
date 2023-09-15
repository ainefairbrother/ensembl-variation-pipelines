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
  example_vep_config = "${projectDir}/assets/example_vep_config.ini"
  
  '''
  if [[ ! -e !{vep_config} || !{force_create_config} == 1 ]]; then
    cp !{example_vep_config} !{vep_config}
    
    #generate_vep_config.py \
    #  !{genome} \
    #  !{version} \
    #  --ini_file !{ini_file} \
    #  --vep_config !{vep_config}
  fi
  '''
}
