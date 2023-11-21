#!/usr/bin/env nextflow

process INDEX_VCF {
  label 'bcftools'
  
  input:
  tuple val(meta), path(vcf)
  
  output:
  tuple val(meta), path(new_vcf), path(vcf_index)
  
  shell:
  index_type = meta.index_type
  flag_index = (index_type == "tbi" ? "-t" : "-c")
  new_vcf = "${meta.genome}-${meta.source}.vcf.gz"
  vcf_index = new_vcf + ".${index_type}"
  
  '''
  ln -sf !{vcf} !{new_vcf}
  bcftools index !{flag_index} !{new_vcf}
  '''
}
