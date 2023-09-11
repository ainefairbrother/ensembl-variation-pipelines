#!/usr/bin/env nextflow

/*
* This script will change the variant id from ClinVar accession to ClinVar variant ids 
*/

process UPDATE_IDS {
  input: 
  tuple path(input_file), val(genome), val(source), val(priority), val(index_type)
  
  output:
  tuple path(input_file), val(genome), val(source), val(priority), val(index_type)
  
  when:
  rename_clinvar_ids = params.rename_clinvar_ids
  output_file = file(input_file).getName().replace("processed", "post_processed")
  
  shell:
  '''
  if [[ !{rename_clinvar_ids} == 1 && !{source} =~ ClinVar ]]; then
    pyenv local variation-eva
    rename_clinvar_ids.py !{input_file} -O !{output_file}
    
    mv !{output_file} !{input_file}
  fi
  '''
}