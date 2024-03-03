#!/usr/bin/env nextflow

process SUMMARY_STATS {
  input: 
  tuple val(meta), path(vcf), path(vcf_index)

  output:
  tuple val(meta), path(output_file), path(vcf_index)
  
  shell:
  output_file =  "UPDATED_SS_" + file(vcf).getName()
  uncomressed_output_file = output_file.replaceFirst(/.gz$/, "")
  index_type = meta.index_type
  flag_index = (index_type == "tbi" ? "-t" : "-c")
  vcf_index = output_file + ".${index_type}"

  '''
  summary_stats.py !{vcf} -O !{uncomressed_output_file}
  bgzip !{uncomressed_output_file}
  bcftools index !{flag_index} !{output_file}
  '''
}