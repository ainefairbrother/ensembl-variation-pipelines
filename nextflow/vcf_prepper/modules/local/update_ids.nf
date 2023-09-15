#!/usr/bin/env nextflow

process UPDATE_IDS {
  input: 
  tuple val(meta), path(vcf)
  
  output:
  tuple val(meta), path(output_file)
  
  when:
  output_file =  "UPDATED_" + file(vcf).getName()
  source = meta.source
  rename_clinvar_ids = params.rename_clinvar_ids
  
  shell:
  '''
  if [[ !{rename_clinvar_ids} == 1 && !{source} =~ ClinVar ]]; then
    pyenv local variation-eva
    rename_clinvar_ids.py !{vcf} -O !{output_file}
  else
    # file must exist as we are using 'path'
    mv !{vcf} !{output_file}
  fi
  '''
}