#!/usr/bin/env nextflow

process SPLIT_VCF_CHR {
  label 'bcftools'
  
  input:
  tuple val(meta), path(vcf), path(index_file), path(chr_file)

  output:
  tuple val(meta), path("split.*.vcf.gz")

  shell:
  
  '''
  chr_file=!{chr_file}
  chr=$(basename ${chr_file/.chrom/})
  
  bcftools view -Oz -o split.${chr}.vcf.gz !{vcf} "${chr}"
  
  rm ${chr_file}
  '''
}
