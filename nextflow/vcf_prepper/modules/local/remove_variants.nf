#!/usr/bin/env nextflow

process REMOVE_VARIANTS {
  label 'bigmem'
  
  input:
  tuple val(meta), path(vcf)
  
  output:
  tuple val(meta), path(output_file)
  
  shell:
  output_file =  "REMOVED_" + file(vcf).getName()
  remove_nonunique_ids = params.remove_nonunique_ids ? "--remove_nonunique_ids" : ""
  remove_patch_regions = params.remove_patch_regions ? "--remove_patch_regions" : ""
  
  '''
  pyenv local variation-eva
  remove_variants.py \
    !{vcf} \
    !{remove_nonunique_ids} \
    !{remove_patch_regions} \
    -O !{output_file}
  '''
}