#!/usr/bin/env nextflow

process READ_VCF_CHR {
  input:
  tuple val(meta), path(vcf), path(vcf_index)

  output:
  tuple val(meta), path(vcf), path(vcf_index), path("*.chrom")

  shell:
  '''
  tabix !{vcf} -l | awk '{print $1".chrom"}' | xargs touch
  '''
}
