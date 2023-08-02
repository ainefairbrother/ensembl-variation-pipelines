#!/usr/bin/env nextflow

/*
* This script concatenates multiple bed files into one
*/

process concatBed {
  input: 
  tuple val(original_vcf), path(bed_files), val(genome), val(source), val(priority)
  
  output: 
  tuple val(original_vcf), path(output_bed), val(genome), val(source), val(priority)
  
  afterScript 'rm all.bed'
  
  shell:
  output_bed = file(original_vcf).getName().replace(".vcf.gz", ".bed")
  
  '''
  cat !{bed_files} > all.bed
  LC_COLLATE=C sort -S1G -k1,1 -k2,2n all.bed > !{output_bed}
  '''
}
