#!/usr/bin/env nextflow

/*
* This script filter out variants from some specific seq region (patches, contigs, etc.)
*/

process filterChr {
  input: 
  tuple val(source), val(vcfFile)
  val genome_dir
  
  output:
  val source, emit: source
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
  output_file_name=$(basename ${vcf_file/_VEP/_filtered_VEP})
  output_file=${output_dir}/${output_file_name}
  
  # seq regions that we want to keep
  keep_regions=$(tabix ${vcf_file} -l | grep -v -E '(CTG|PATCH|TEST)' |  paste -sd,)
  
  bcftools view -r ${keep_regions} ${vcf_file} -Oz -o ${output_file}
  '''
}