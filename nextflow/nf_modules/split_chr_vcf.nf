#!/usr/bin/env nextflow

/* 
 * Script to split a VCF file into multiple smaller VCFs
 */

nextflow.enable.dsl=2

// defaults
params.cpus = 1

process splitChrVCF {
  input:
  tuple path(vcf), path(index_file), path(chr_file), val(genome)

  output:
  tuple val("${vcf}"), path("split.*.vcf.gz"), val(genome)
  
  afterScript 'rm ${chr_file}'

  shell:
  '''
  chr_file=!{chr_file}
  chr=$(basename ${chr_file/.chrom/})
  
  bcftools view -Oz -o split.${chr}.vcf.gz !{vcf} "${chr}"
  '''
}
