#!/usr/bin/env nextflow

/*
* This script generate bed file from VCF file
*/

process vcfToBed {
  input: 
  val config_done
  tuple val(original_vcf), path(vcf_file), val(genome), val(source), val(priority)
  
  output:
  tuple val(original_vcf), path(output_filename), val(genome), val(source), val(priority)
  
  shell:
  output_filename = file(original_vcf).getSimpleName() + "-" + vcf_file
  rank_file = params.rank_file
  
  '''
  !{projectDir}/../../bin/vcf_to_bed \
    !{vcf_file} \
    !{output_filename} \
    !{rank_file}
    
  rm !{vcf_file}
  '''
}
