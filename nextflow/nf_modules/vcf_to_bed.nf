#!/usr/bin/env nextflow

/*
* This script generate bed file from VCF file
*/

process vcfToBed {
  input: 
  tuple val(original), path(vcfFile), path(indexFile)
  
  output:
  tuple val(original), path("${original}-${vcfFile}.bed"), emit: bed
  
  shell:
  '''
  output_filename=!{original}-!{vcfFile}.bed
  
  # run the track generation script
  !{projectDir}/../../src/rust/ensembl/vcf_to_bed/bin/vcf_to_bed \
    !{vcfFile} \
    ${output_filename} \
    !{projectDir}/../nf_config/variation_consequnce_rank.json
    
  rm !{vcfFile}
  '''
}
