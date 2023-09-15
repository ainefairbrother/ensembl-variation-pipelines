#!/usr/bin/env nextflow

process CONCAT_BEDS {
  input: 
  tuple val(meta), path(bed_files)
  
  output: 
  tuple val(meta), path(output_bed)
  
  afterScript 'rm all.bed'
  
  shell:
  source = meta.source
  output_bed =  "variant-${source}.bed"
    
  '''
  cat !{bed_files} > all.bed
  LC_COLLATE=C sort -S1G -k1,1 -k2,2n all.bed > !{output_bed}
  '''
}
