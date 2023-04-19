#!/usr/bin/env nextflow

/*
* This script filter out variants from some specific seq region (patches, contigs, etc.)
*/

process indexVCF {
  input:
  val input_vcf
  
  output:
  tuple val(input_vcf), val("${input_vcf}.gz")
  
  shell:
  '''
  # indexing the file
  bcftools index -t !{input_vcf}
  '''
}
