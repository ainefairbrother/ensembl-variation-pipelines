#!/usr/bin/env nextflow

process RENAME_CHR {
  label 'bcftools'
  
  input:
  tuple val(meta), val(vcf), val(vcf_index)
  
  output:
  tuple val(meta), path(output_file)
  
  shell:
  output_file =  "RENAMED_" + file(vcf).getName()
  synonym_file = meta.synonym_file
  
  '''
  bcftools annotate --no-version --force --rename-chrs !{synonym_file} !{vcf} -Oz -o !{output_file}
  '''
}
