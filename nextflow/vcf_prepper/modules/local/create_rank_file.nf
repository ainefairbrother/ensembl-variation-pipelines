#!/usr/bin/env nextflow

process CREATE_RANK_FILE {
  input:
  val rank_file
  
  output: 
  val rank_file
  
  shell:
  '''
  if [[ -f !{rank_file} ]]; then
    rm !{rank_file}
  fi
  
  generate_consequence_rank.pl -o !{rank_file}
  '''
}
