#!/usr/bin/env nextflow

/*
* This script filter out variants from some specific seq region (patches, contigs, etc.)
*/

process filterChr {
  input:
  val renamedVCFFile
  
  output:
  env output_file, emit: vcfFile
  
  shell:
  '''
  renamed_vcf_file=!{renamedVCFFile}
  
  # get output dir
  output_dir=$(dirname ${renamed_vcf_file})
  
  # format output file
  output_file_name=$(basename ${renamed_vcf_file/_renamed_VEP/_filtered_VEP})
  output_file_name=$(basename ${output_file_name/.gz/})    # we won't compress the output file
  output_file=${output_dir}/${output_file_name}
  
  # seq regions that we want to keep
  keep_regions=$(tabix ${renamed_vcf_file} -l | grep -v -E '(CTG|PATCH|TEST)' |  paste -sd,)
  
  # only keep the regions that we want to keep
  bcftools view -r ${keep_regions} ${renamed_vcf_file} -Oz -o ${output_file}
  
  rm ${renamed_vcf_file} ${renamed_vcf_file}.tbi
  '''
}
