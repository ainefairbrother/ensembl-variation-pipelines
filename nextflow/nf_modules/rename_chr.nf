#!/usr/bin/env nextflow

/*
* This script rename the seq region synonyms in the VCF file
*/

process renameChr {
  input:
  val input_vcf
  path synonym_file
  
  output:
  tuple env(output_dir), env(prefix)
  
  shell:
  '''
  # format input and output file name
  input_file=!{input_vcf}
  output_file=${input_file/.vcf.gz_vep_out.vcf.gz/_renamed_VEP.vcf.gz}
  
  # rename chr synonyms
  bcftools annotate --no-version --force --rename-chrs !{synonym_file} ${input_file} -Oz -o ${output_file}
  
  # for next job create parameters
  output_dir=$(dirname ${output_file})
  prefix=$(basename ${output_file/_renamed_VEP.vcf.gz/})
  '''
}
