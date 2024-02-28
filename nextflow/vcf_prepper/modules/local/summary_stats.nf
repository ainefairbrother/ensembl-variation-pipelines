#!/usr/bin/env nextflow

process SUMMARY_STATS {
  input: 
  tuple val(meta), path(vcf), path(vcf_index)

  output:
  tuple val(meta), path(output_file), path(vcf_index)
  
  shell:
  output_file =  "UPDATED_" + file(vcf).getName()
  index_type = meta.index_type
  flag_index = (index_type == "tbi" ? "-t" : "-c")
  vcf_index = output_file + ".${index_type}"

  '''
  summary_stats.py !{vcf} -O !{output_file}
  bcftools index !{flag_index} !{output_file}
  '''
}