#!/usr/bin/env nextflow

/*
* This script create config files required for later stage
*/

process CREATE_CONFIGS {
  input:
  val rank_file_done
  val genome
  
  output: 
  val "done"
  
  shell:
  ini_file = params.ini_file
  version = params.version
  synonym_file = "${moduleDir}/../../common_files/${genome}.txt"
  chrom_sizes = "${moduleDir}/../../common_files/${genome}.chrom.sizes"
  skip_create_config = params.skip_create_config
  
  '''
  if [[ !{skip_create_config} == 0 ]]; then
    generate_configs.py \
      !{genome} \
      !{version} \
      --ini_file !{ini_file} \
      --synonym_file !{synonym_file} \
      --chrom_sizes !{chrom_sizes}
  fi
  '''
}
