#!/usr/bin/env nextflow

/*
* This script create config files required for later stage
*/

process createConfigs {
  input:
  val rank_file_done
  val genome
  
  output: 
  val "done"
  
  shell:
  ini_file = params.ini_file
  version = params.version
  skip_create_config = params.skip_create_config
  
  '''
  if [[ !{skip_create_config} == 0 ]]; then
    python !{moduleDir}/../../src/python/ensembl/scripts/generate_configs.py \
      -I !{ini_file}  \
      --genome !{genome} \
      --version !{version}
  fi
  '''
}
