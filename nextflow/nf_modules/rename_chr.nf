#!/usr/bin/env nextflow

/*
* This script rename the seq region synonyms in the VCF file
*/

process renameChr {
  input:
  tuple val(source), val(vcfFile)
  val genome_dir
  
  output:
  env output_file, emit: vcfFile
  
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
  
  # format output file name
  output_file_name=$(basename ${vcf_file/_VEP/_renamed_VEP})
  output_file=${output_dir}/${output_file_name}
  
  # synonym file that contains mapping between seq region name and synonym
  chr_synonym_file=!{projectDir}/../nf_config/homo_sapiens_grch38_synonyms.txt
  
  # rename chr synonyms
  bcftools annotate --no-version --rename-chrs ${chr_synonym_file} ${vcf_file} -Oz -o ${output_file}
  
  # indexing the file because next stage requires it
  bcftools index -t ${output_file}
  '''
}
