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
    output_filename=$(basename ${vcf_file_template/\\#\\#CHR\\#\\#/all})
    temp_output_filename=temp_${output_filename/.gz/}
    
    if [[ "${source}" == "10000_Genomes" ]]; then
      for file in $input_files;
      do
        # remove genotype from name as we are actually removing genotypes
        new_file=$(basename ${file/genotypes./})
        new_file=${new_file/.gz/}

        # get the headers removing all genotype related content
        tabix ${file} -H | grep -v "^##FORMAT" | awk '{if(/^#CHR/) for(i=1;i<=8;i++) printf $i"\t"; else print}' > ${output_dir}/filtered_${new_file}
        
        # get the content of the file
        chr=$(tabix ${file} -l)
        echo "" >> ${output_dir}/filtered_${new_file}
        tabix ${file} ${chr} | cut -d'	' -f1-8 >> ${output_dir}/filtered_${new_file}
        
        # compress and tabix
        bgzip ${output_dir}/filtered_${new_file}
        tabix -p vcf ${output_dir}/filtered_${new_file}.gz
      done
      
      input_files=${output_dir}/filtered_${new_file/chr[1-9XY]/*}.gz
      
      # vcf-concat cannot parse bgzip file - use bcftools and hope it works 
      # vcf-concat --pad-missing ${input_files} > ${output_dir}/${temp_output_filename}
    fi
      
    # concat the input files
    bcftools concat --no-version ${input_files} -Oz -o ${output_dir}/${temp_output_filename}
    
    # sort the output file and compress it
    temp_dir=${output_dir}/temp_${source}
    mkdir -p ${temp_dir}
    bcftools sort -T ${temp_dir} ${output_dir}/${temp_output_filename} -Oz -o ${output_dir}/${output_filename}
    
    # delete the temporary dirs
    rm ${output_dir}/${temp_output_filename}
    
    # delete the temporary input files if the source is 1000 Genomes
    if [[ "${source}" == "10000_Genomes" ]]; then
      # not using the $input_files directly as in weird case some actual input source gets deleted
      rm ${output_dir}/filtered_*
    fi
  fi
  '''
}