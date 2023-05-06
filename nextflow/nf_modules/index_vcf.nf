#!/usr/bin/env nextflow

/*
* This script filter out variants from some specific seq region (patches, contigs, etc.)
*/

process indexVCF {
  input:
  tuple val(input_vcf), val(genome), val(source), val(priority), val(index_type)
  
  output:
  tuple path(input_vcf), path("${input_vcf}.{tbi,csi}"), val(genome), val(source), val(priority)
  
  shell:
  '''
  if [[ "!{index_type}" == "tbi" ]]; then
    bcftools index -t !{input_vcf}
  else
    bcftools index -c !{input_vcf}
  fi
  '''
}
