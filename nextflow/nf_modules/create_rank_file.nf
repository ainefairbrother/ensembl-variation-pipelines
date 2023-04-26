#!/usr/bin/env nextflow

/*
* This script create variation consequence rank file in json format using the Ensembl Variation API
*/

process createRankFile {
  input:
  val rank_file
  output: 
  val "done"
  
  shell:
  '''
  if [[ -f !{rank_file} ]]; then
    rm !{rank_file}
  fi
  
  perl !{projectDir}/../../src/perl/ensembl/scripts/generate_consequence_rank.pl \
    -o !{rank_file}
  '''
}
