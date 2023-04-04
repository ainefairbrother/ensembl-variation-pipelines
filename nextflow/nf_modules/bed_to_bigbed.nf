#!/usr/bin/env nextflow

/*
* This script generate bed file from VCF file
*/

process bedToBigBed {
  input: 
  tuple val(original), path(bedFiles)
  
  shell:
  '''
  cat !{bedFiles} > all.bed
  #LC_COLLATE=C sort -S1G -k1,1 -k2,2n all.bed > all.sorted.bed
  #bedToBigBed -type=bed3+6 all.sorted.bed !{projectDir}/../nf_config/chrom_sizes/homo_sapiens_grch38.chrom.sizes all.bb
  
  #rm all.bed
  '''
}
