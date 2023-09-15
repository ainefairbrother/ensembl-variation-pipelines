#!/usr/bin/env nextflow

process INDEX_VCF {
  label 'bcftools'
  
  input:
  tuple val(meta), val(vcf)
  
  output:
  tuple val(meta), val(vcf), val(vcf_index)
  
  shell:
  index_type = meta.index_type
  flag_index = (index_type == "tbi" ? "-t" : "-c")
  vcf_index = file(vcf) + ".${index_type}"
  
  '''
  bcftools index !{flag_index} !{vcf}
  '''
}
