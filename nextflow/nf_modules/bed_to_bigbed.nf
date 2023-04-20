#!/usr/bin/env nextflow

/*
* This script generate bed file from VCF file
*/

process bedToBigBed {
  input: 
  tuple val(original), path(bed_files), val(genome)
  
  afterScript 'rm all.bed'
  
  shell:
  '''
  chrom_sizes=!{projectDir}/../nf_config/chrom_sizes/!{genome}.chrom.sizes
  
  cat !{bed_files} > all.bed
  LC_COLLATE=C sort -S1G -k1,1 -k2,2n all.bed > all.sorted.bed
  bedToBigBed -type=bed3+6 all.sorted.bed ${chrom_sizes} all.bb
  
  
  '''
}
