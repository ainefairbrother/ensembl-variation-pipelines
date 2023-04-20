#!/usr/bin/env nextflow

/*
* This script rename the seq region synonyms in the VCF file
*/

process renameChr {
  input:
  val input_file
  val genome
  
  output:
  tuple env(output_dir), env(prefix), val(genome), env(index_type)
  
  shell:
  '''
  # format input and output file name
  input_file=!{input_file}
  output_file=${input_file/_VEP.vcf.gz/_renamed_VEP.vcf.gz}
  synonym_file=!{projectDir}/../nf_config/synonyms/!{genome}.txt
  
  # rename chr synonyms
  bcftools annotate --no-version --force --rename-chrs ${synonym_file} ${input_file} -Oz -o ${output_file}
  
  # for next job create parameters
  output_dir=$(dirname ${output_file})
  prefix=$(basename ${output_file/_renamed_VEP.vcf.gz/})
  if [[ -f ${input_file}.tbi ]]; then
    index_type=tbi
  else
    index_type=csi
  fi
  '''
}
