#!/usr/bin/env nextflow

/*
* This script post prcess a VCF file after running through VEP. It will -
* - remove any seq region with PATCH/TEST/CTG in it's name
* - remove any variant with duplicated rsIDs
*/

process removeDupIDs {
  input: 
  tuple val(output_dir), val(prefix)
  
  output:
  env output_file
  
  shell: 
  '''
  # format input and output file name
  input_file=!{output_dir}/!{prefix}_renamed_VEP.vcf.gz
  output_file=!{output_dir}/!{prefix}_processed_VEP.vcf.gz
  
  pyenv local variation-eva
  python3 !{projectDir}/../../src/python/ensembl/scripts/remove_duplicate_ids.py ${input_file}
  '''
}