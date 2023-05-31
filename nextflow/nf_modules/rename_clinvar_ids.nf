#!/usr/bin/env nextflow

/*
* This script will change the variant id from ClinVar accession to ClinVar variant ids 
*/

process renameClinvarIDs {
  input: 
  tuple val(input_file), val(genome), val(source), val(priority), val(index_type)
  
  output:
  tuple val(input_file), val(genome), val(source), val(priority), val(index_type)
  
  when:
  rename_clinvar_ids = params.rename_clinvar_ids
  
  shell:
  '''
  if [[ !{rename_clinvar_ids} == 1 && source =~ ClinVar ]]; then
    # format input and output file name
    input_file=!{input_file}
    output_file=${input_file/_processed_VEP.vcf.gz/_post_processed_VEP.vcf.gz}
    
    pyenv local variation-eva
    python3 !{moduleDir}/../../src/python/ensembl/scripts/rename_clinvar_ids.py ${input_file} ${output_file}
    
    mv ${output_file} ${input_file}
  fi
  '''
}