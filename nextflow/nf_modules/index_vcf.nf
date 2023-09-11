#!/usr/bin/env nextflow

/*
* This script filter out variants from some specific seq region (patches, contigs, etc.)
*/

process indexVCF {
  label 'bcftools'
  
  input:
  tuple val(input_vcf), val(genome), val(source), val(priority), val(index_type)
  
  output:
  tuple val(input_vcf), env(vcf_index), val(genome), val(source), val(priority)
  
  shell:
  '''
  if [[ "!{index_type}" == "tbi" ]]; then
    bcftools index -t !{input_vcf}
    vcf_index=!{input_vcf}.tbi
  else
    bcftools index -c !{input_vcf}
    vcf_index=!{input_vcf}.csi
  fi
  '''
}
