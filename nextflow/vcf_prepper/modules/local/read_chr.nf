#!/usr/bin/env nextflow

/* 
 * Script to split a VCF file into multiple smaller VCFs
 */

process READ_CHR {
  input:
  tuple path(vcf), path(vcf_index), val(genome), val(source), val(priority)

  output:
  tuple path(vcf), path(vcf_index), path("*.chrom"), val(genome), val(source), val(priority)

  shell:
  '''
  tabix !{vcf} -l | awk '{print $1".chrom"}' | xargs touch
  '''
}
