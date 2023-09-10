#!/usr/bin/env nextflow

/*
* This script rename the seq region synonyms in the VCF file
*/

process RENAME_CHR {
  label 'bcftools'
  
  input:
  val input_file
  val genome
  val source
  val priorities
  
  output:
  tuple path(output_file), val(genome), val(source), val(priority), env(index_type)
  
  shell:
  priority = priorities[genome][source]
  output_file = file(input_file).getName().replace("_VEP.vcf.gz", "_renamed_VEP.vcf.gz")
  synonym_file = "${moduleDir}/../../common_files/synonyms/${genome}.txt"
  
  '''
  # rename chr synonyms
  bcftools annotate --no-version --force --rename-chrs !{synonym_file} !{input_file} -Oz -o !{output_file}
  
  # for next job create parameters
  if [[ -f !{input_file}.tbi ]]; then
    index_type=tbi
  else
    index_type=csi
  fi
  '''
}
