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
  temp_dir = "tmp"
    
  '''
  # tmp directory for sorting big files
  mkdir -p !{temp_dir}

  cat !{bed_files} > all.bed
  LC_COLLATE=C sort -T !{temp_dir} -S1G -k1,1 -k2,2n all.bed > !{output_bed}

  rm -r !{temp_dir}
  '''
}
