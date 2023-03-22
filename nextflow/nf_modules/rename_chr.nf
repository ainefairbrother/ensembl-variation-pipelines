#!/usr/bin/env nextflow

/*
* This script rename the seq region synonyms in the VCF file
*/

process renameChr {
  input: 
  val source
  val vcfFile
  val genome_dir
  
  shell:
  '''
  vcf_file=!{vcfFile}
  
  # remove any whitespace in source name
  source="!{source}"
  source=${source// /_}
  
  # create a source dir under output dir
  source_dir=!{genome_dir}/${source}
  mkdir -p $source_dir
  
  # create output dir
  output_dir=${source_dir}/vcfs
  mkdir -p $output_dir
  
  # synonym file that contains mapping between seq region name and synonym
  chr_synonym_file=!{projectDir}/../nf_config/homo_sapiens_grch38_synonyms.txt
  
  bcftools annotate --rename-chrs ${chr_synonym_file} ${vcf_file} -Oz -o ${vcf_file}
  '''
}