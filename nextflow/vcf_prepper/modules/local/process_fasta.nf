#!/usr/bin/env nextflow

process PROCESS_FASTA {
  cache false
  
  input:
  val meta

  output:
  val genome
  
  shell:
  genome = meta.genome
  species = meta.species
  assembly = meta.assembly
  version = params.version
  ini_file = params.ini_file
  fasta_dir = meta.fasta_dir
  
  '''
  process_fasta.py \
    !{species} \
    !{assembly} \
    !{version} \
    --ini_file !{ini_file} \
    --fasta_dir !{fasta_dir}
  '''
}
