#!/usr/bin/env nextflow

/*
* This script post prcess a VCF file after running through VEP. It will -
* - remove any seq region with PATCH/TEST/CTG in it's name
* - remove any variant with duplicated rsIDs
*/

process REMOVE_VARIANTS {
  label 'bigmem'
  
  input: 
  tuple path(input_file), val(genome), val(source), val(priority), val(index_type)
  
  output:
  tuple path(output_file), val(genome), val(source), val(priority), val(index_type)
  
  shell:
  remove_nonunique_ids = params.remove_nonunique_ids ? "--remove_nonunique_ids" : ""
  remove_patch_regions = params.remove_patch_regions ? "--remove_patch_regions" : ""
  output_file = file(input_file).getName().replace("renamed", "processed")
  
  '''
  pyenv local variation-eva
  remove_variants.py \
    !{input_file} \
    !{remove_nonunique_ids} \
    !{remove_patch_regions} \
    -O !{output_file}
  '''
}