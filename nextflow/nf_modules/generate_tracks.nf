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
  
  file_type=$(file ${vcf_file})
  if [[ ${file_type} =~ compressed ]]; then
    temp_dir=${PWD}/${source// /_}_temp
    mkdir -p ${temp_dir}
    
    vcf_filename=$(basename ${vcf_file})
    uncompressed_vcf=${temp_dir}/${vcf_filename/.gz/}
    bgzip -c -d ${vcf_file} > ${uncompressed_vcf}
    
    vcf_file=${uncompressed_vcf}
  fi
  
  out_dir=${PWD}/${source}_out
  mkdir -p ${out_dir}
  
  config_dir=!{projectDir}/../nf_config
  
  perl !{projectDir}/../../src/perl/ensembl/scripts/generate_tracks.pl ${vcf_file} ${out_dir} ${config_dir}
  
  if [[ -d ${temp_dir} ]]; then
    rm -r ${temp_dir}
  fi
  '''
}
