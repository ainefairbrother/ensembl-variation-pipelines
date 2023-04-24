#!/usr/bin/env nextflow

/* 
 * Script to split a VCF file into multiple smaller VCFs
 */

nextflow.enable.dsl=2

// defaults
params.cpus = 1

process readChrVCF {
  input:
  tuple path(vcf), path(vcf_index), val(genome), val(source)

  output:
  tuple path(vcf), path(vcf_index), path("*.chrom"), val(genome), val(source)

  shell:
  '''
  tabix !{vcf} -l | awk '{print $1".chrom"}' | xargs touch
  '''
}
