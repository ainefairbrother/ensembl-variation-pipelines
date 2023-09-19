#!/usr/bin/env nextflow

process GENERATE_CHROM_SIZES {
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
  chrom_sizes = meta.chrom_sizes
  
  '''
  # for -resume ability we are creating the files in the work dir instead of giving --chrom_sizes param
  if [[ ! -e !{chrom_sizes} || !{force_create_config} == 1 ]]; then
    generate_chrom_sizes.py \
      !{species} \
      !{assembly} \
      !{version} \
      --ini_file !{ini_file} \
      --chrom_sizes !{chrom_sizes}
  fi
  '''
}
