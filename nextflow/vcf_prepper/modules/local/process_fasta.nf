#!/usr/bin/env nextflow

process PROCESS_FASTA {
  cache false
  
  input:
  val meta

  output:
  val genome
  
  shell:
  genome = meta.genome
  version = params.version
  ini_file = params.ini_file
  fasta_dir = params.fasta_dir
  
  '''
  process_fasta.py \
    !{genome} \
    !{version} \
    --ini_file !{ini_file} \
    --fasta_dir !{fasta_dir}
  '''
}
