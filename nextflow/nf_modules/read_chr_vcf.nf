#!/usr/bin/env nextflow

/* 
 * Script to split a VCF file into multiple smaller VCFs
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"
params.outdir = ""
params.cpus = 1

process readChrVCF {
  input:
  path(vcf)
  path(indexFile)

  output:
  tuple path(vcf), path(indexFile), path("*.chrom")

  shell:
  '''
  tabix !{vcf} -l | awk '{print $1".chrom"}' | xargs touch
  '''
}
