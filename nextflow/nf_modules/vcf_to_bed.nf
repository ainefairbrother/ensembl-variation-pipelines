#!/usr/bin/env nextflow

/*
* This script generate bed file from VCF file
*/

process vcfToBed {
  input: 
  tuple val(original_vcf), path(vcf_file), val(genome)
  
  output:
  tuple val(original_vcf), path(output_filename), val(genome)
  
  shell:
  output_filename = file(original_vcf).getSimpleName() + "-" + vcf_file
  '''
  # run the track generation script
  !{projectDir}/../../bin/vcf_to_bed \
    !{vcf_file} \
    !{output_filename} \
    !{projectDir}/../nf_config/variation_consequnce_rank.json
    
  rm !{vcf_file}
  '''
}
