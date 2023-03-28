#!/usr/bin/env nextflow

/*
* This script create combined focus VCF file from multiple sources
*/

process createFocusVCF {
  input:
  path vcfFiles
  path indexFiles
  val genome_dir
  
  shell:
  '''
  # create a source dir under output dir
  source_dir=!{genome_dir}/focus
  mkdir -p $source_dir
  
  # create output dir
  output_dir=${source_dir}/vcfs
  mkdir -p $output_dir
  
  # format output file name
  output_file_name=$(basename !{genome_dir}).vcf.gz
  output_file=${output_dir}/${output_file_name}
  
  # concat all the vcf vcf files
  bcftools concat --no-version -a -D !{vcfFiles} -Oz -o ${output_file}_tmp
  
  # sort the file
  mkdir -p temp
  bcftools sort -T temp ${output_file}_tmp -Oz -o ${output_file}
  rm ${output_file}_tmp
  
  # index the file
  bcftools index -t ${output_file}
  '''
}