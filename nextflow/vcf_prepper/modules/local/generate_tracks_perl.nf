#!/usr/bin/env nextflow

/*
 * See the NOTICE file distributed with this work for additional information
 * regarding copyright ownership.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
 
process GENERATE_TRACKS_PERL {
  input: 
  tuple val(source), val(vcfFile)
  
  shell:
  '''
  vcf_file="!{vcfFile}"
  source="!{source}"
  
  # remove any whitespace in source name
  source=${source// /_}
  
  # if vcf file is uncompressed we need to compress it
  temp_dir=${PWD}/${source}_temp
  mkdir -p ${temp_dir}  
  file_type=$(file ${vcf_file})
  if [[ ! ${file_type} =~ compressed ]]; then
    # compress vcf file
    vcf_filename=$(basename ${vcf_file})
    compressed_vcf=${temp_dir}/${vcf_filename/.gz/}
    bgzip -c ${vcf_file} > ${compressed_vcf}
    
    # set the vcf_file to compressed vcf file
    vcf_file=${compressed_vcf}
  fi
  
  # tabix index the compressed vcf file if not already exist
  tabix ${vcf_file} -l > /dev/null 2>&1 || tabix -p vcf -f ${vcf_file}
  
  out_dir=${PWD}/${source}_out
  mkdir -p ${out_dir}
  
  config_dir=!{projectDir}/assets
  
  # run the track generation script
  generate_tracks.pl ${vcf_file} ${out_dir} ${config_dir}

  # remove temp dir with compressed vcf file
  if [[ -d ${temp_dir} ]]; then
    rm -r ${temp_dir}
  fi
  '''
}
