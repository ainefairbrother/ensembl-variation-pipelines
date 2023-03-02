#!/usr/bin/env nextflow

/*
* This script generate bigBed and bigwig files from VCF file
*/

process generateTracks {
  input: 
  tuple val(source),val(vcfFile)
  
  shell:
  '''
  vcf_file="!{vcfFile}"
  source="!{source}"
  
  # remove any whitespace in source name
  source=${source// /_}
  
  # if vcf file is compressed we need to uncompress it
  file_type=$(file ${vcf_file})
  if [[ ${file_type} =~ compressed ]]; then
    # create a temp dir to put uncompressed vcf file - we will remove it later
    temp_dir=${PWD}/${source}_temp
    mkdir -p ${temp_dir}
  
    # uncompress vcf file
    vcf_filename=$(basename ${vcf_file})
    uncompressed_vcf=${temp_dir}/${vcf_filename/.gz/}
    bgzip -c -d ${vcf_file} > ${uncompressed_vcf}
    
    # set the vcf_file to uncompressed vcf file
    vcf_file=${uncompressed_vcf}
  fi
  
  out_dir=${PWD}/${source}_out
  mkdir -p ${out_dir}
  
  config_dir=!{projectDir}/../nf_config
  
  # run the track generation script
  perl !{projectDir}/../../src/perl/ensembl/scripts/generate_tracks.pl ${vcf_file} ${out_dir} ${config_dir}

  # remove temp dir with uncompressed vcf file
  if [[ -d ${temp_dir} ]]; then
    rm -r ${temp_dir}
  fi
  '''
}
