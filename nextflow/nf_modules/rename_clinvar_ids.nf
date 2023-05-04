#!/usr/bin/env nextflow

/*
* This script will change the variant id from ClinVar accession to ClinVar variant ids 
*/

process renameClinvarIDs {
  input: 
  tuple val(output_dir), val(prefix), val(genome), val(source), val(priority), val(index_type)
  
  output:
  tuple env(input_file), val(prefix), val(genome), val(source), val(priority), val(index_type)
  
  when:
  params.rename_clinvar_ids == 1 && source =~ /ClinVar/
  
  shell:
  rename_clinvar_ids = params.rename_clinvar_ids
  
  '''
  # format input and output file name
  input_file=!{output_dir}/!{prefix}_processed_VEP.vcf.gz
  output_file=!{output_dir}/!{prefix}_processed_bkp_VEP.vcf.gz
  
  pyenv local variation-eva
  python3 !{projectDir}/../../src/python/ensembl/scripts/rename_clinvar_ids.py ${input_file} !{output_file}
  
  mv !{output_file} !{input_file}
  '''
}