#!/usr/bin/env nextflow

/* 
 * Script to split a VCF file into multiple smaller VCFs
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"
params.outdir = ""
params.cpus = 1

process splitChrVCF {
  input:
  tuple path(vcf), path(indexFile), path(chrFile)

  output:
  tuple val("${vcf}"), path("split.*.vcf.gz"), path(indexFile), emit: files

  shell:
  '''
  chr_file=!{chrFile}
  chr=$(basename ${chr_file/.chrom/})
  bcftools view -O z -o split.${chr}.vcf.gz !{vcf} "${chr}"
  '''
}
