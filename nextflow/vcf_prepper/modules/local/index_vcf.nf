#!/usr/bin/env nextflow

/*
* This script filter out variants from some specific seq region (patches, contigs, etc.)
*/

process INDEX_VCF {
  label 'bcftools'
  
  input:
  tuple path(input_vcf), val(genome), val(source), val(priority), val(index_type)
  
  output:
  tuple path(input_vcf), path(vcf_index), val(genome), val(source), val(priority)
  
  shell:
  vcf_index = file(input_vcf).getName() + ".tbi"
  
  '''
  if [[ "!{index_type}" == "tbi" ]]; then
    bcftools index -t !{input_vcf}
  else
    bcftools index -c !{input_vcf}
    !{vcf_index}=!{input_vcf}.csi
  fi
  '''
}
