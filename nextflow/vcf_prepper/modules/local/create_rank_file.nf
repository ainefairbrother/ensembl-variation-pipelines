#!/usr/bin/env nextflow

/*
* This script create variation consequence rank file in json format using the Ensembl Variation API
*/

process CREATE_RANK_FILE {
  input:
  val rank_file
  
  output: 
  val "done"
  
  shell:
  '''
  if [[ -f !{rank_file} ]]; then
    rm !{rank_file}
  fi
  
  generate_consequence_rank.pl -o !{rank_file}
  '''
}
