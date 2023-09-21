#!/usr/bin/env nextflow

process DOWNLOAD_SOURCE {
  cache false
  
  input:
  tuple val(meta), val(vcf)

  output:
  tuple val(meta), val(output_vcf)
  
  shell:
  file_type = meta.file_type
  output_vcf = file_type == "remote" ? meta.genome_temp_dir + "/" + file(vcf).getName() : vcf
  
  '''
  if [[ !{file_type} == "remote" ]]; then
    wget !{vcf} -O !{output_vcf}
  fi
  '''
}
