#!/usr/bin/env nextflow

/*
* This script rename all the duplicated rsID to '.' in a VCF file 
*/

process removeDupIDs {
  input:
  val filteredVCFFile
  
  output:
  env output_file, emit: vcfFile
  env index_file, emit: indexFile
  
  shell:
  '''
  filtered_vcf_file=!{filteredVCFFile}
  
  # get output dir
  output_dir=$(dirname ${filtered_vcf_file})
  
  # format output file name
  output_file_name=$(basename ${filtered_vcf_file/_filtered_VEP/_processed_VEP})
  output_file=${output_dir}/${output_file_name}
  
  # we will do infile change - keep filtered_vcf_file as backup in case the job fails at mid-point
  cp ${filtered_vcf_file} ${output_file}
  
  # get the duplicated rsID list
  bcftools view --no-version ${output_file} | awk '!/^#/{if (uniq[$3]++) print($3)}' | sort -u > duplicated_ids.txt 
  
  # only keep the regions that we want to keep
  while IFS= read -r id;
  do
    awk -i inplace -v var=$id -F'\t' 'BEGIN {OFS = FS} {gsub(var, ".", $3)}1' ${output_file}
  done < duplicated_ids.txt
  
  # bgzip and sort for next step
  bgzip ${output_file}
  bcftools index -t ${output_file}.gz
  
  #output_file=${output_file}.gz
  index_file=${output_file}.tbi
  
  rm ${filtered_vcf_file}
  '''
}