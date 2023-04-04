#!/usr/bin/env nextflow

/*
* This script rename the seq region synonyms in the VCF file
*/

process renameChr {
  input:
  tuple val(original), path(vcfFile), path(indexFile)
  
  output:
  tuple val(original), path("renamed-${original}-*.vcf.gz"), path("renamed-${original}-*.vcf.gz.tbi"), emit: file
  
  shell:
  '''
  # format output file name
  output_file=renamed-!{original}-!{vcfFile}
  
  # synonym file that contains mapping between seq region name and synonym
  chr_synonym_file=!{projectDir}/../nf_config/homo_sapiens_grch38_synonyms.txt
  
  # rename chr synonyms
  bcftools annotate --no-version --force --rename-chrs ${chr_synonym_file} !{vcfFile} -Oz -o ${output_file}
  
  # indexing the file because next stage requires it
  bcftools index -t ${output_file}
  '''
}
