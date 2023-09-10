#!/usr/bin/env nextflow

/* 
 * Script to split a VCF file into multiple smaller VCFs
 */

process SPLIT_VCF {
  label 'bcftools'
  
  input:
  tuple path(vcf), path(index_file), path(chr_file), val(genome), val(source), val(priority)

  output:
  tuple val("${vcf}"), path("split.*.vcf.gz"), val(genome), val(source), val(priority)

  shell:
  '''
  chr_file=!{chr_file}
  chr=$(basename ${chr_file/.chrom/})
  
  bcftools view -Oz -o split.${chr}.vcf.gz !{vcf} "${chr}"
  
  rm ${chr_file}
  '''
}
