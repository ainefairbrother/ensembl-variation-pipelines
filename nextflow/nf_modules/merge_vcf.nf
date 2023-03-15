#!/usr/bin/env nextflow

/*
* This script combines chromosome separated vcf files together
*/

process mergeVCF {
  input: 
  tuple val(source), val(vcfFile)
  val output_dir
  
  shell:
  '''
  vcf_file_template="!{vcfFile}"
  
  # only run if the template have ##CHR## in the name (we are escaping # twice - once from groovy and another for bash)
  if [[ ${vcf_file_template} =~ \\#\\#CHR\\#\\# ]]; then 
    # remove any whitespace in source name
    source="!{source}"
    source=${source// /_}
    
    # input file with wildcard
    input_files=${vcf_file_template/\\#\\#CHR\\#\\#/*}
    
    # create a source dir under output dir
    output_dir=!{output_dir}/${source}
    mkdir -p $output_dir
    
    # configure the output file names
    output_filename=$(basename ${vcf_file_template/\\#\\#CHR\\#\\#/})
    temp_output_filename=temp_${output_filename/.gz/}
    
    if [[ ${vcf_file_template} =~ bgz$ ]]; then
      # vcf-concat cannot parse bgzip file - use bcftools and hope it works 
      bcftools concat --no-version ${input_files} -Oz -o ${output_dir}/${temp_output_filename}
    else
      # we are using vcf-concat because some sources may need padding for missing column (Y chr for 1000G)
      vcf-concat --pad-missing ${input_files} > ${output_dir}/${temp_output_filename}
    fi
    
    # sort the output file and compress it (we create a temp dir for this)
    temp_dir=${output_dir}/temp_${source}
    mkdir -p ${temp_dir}
    # we still have header problem from vcf-concat - it does not include all header lines from all files and some header lines may get missing. bcftools sort 
    # complains if it does not find a header value for a field
    bcftools sort -T ${temp_dir} ${output_dir}/${temp_output_filename} -Oz -o ${output_dir}/${output_filename}
    
    # delete the temporary dirs
    rm $temp_dir
    rm ${output_dir}/${temp_output_filename}
  fi
  '''
}