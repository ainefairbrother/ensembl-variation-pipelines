#!/usr/bin/env nextflow

process GENERATE_CHROM_SIZES {
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
  chrom_sizes = meta.chrom_sizes
  force_create_config = params.force_create_config ? "--force" : ""
  
  '''
  generate_chrom_sizes.py \
    !{species} \
    !{assembly} \
    !{version} \
    --ini_file !{ini_file} \
    --chrom_sizes !{chrom_sizes} \
    !{force_create_config}
  '''
}
