#!/usr/bin/env nextflow

/*
* This script generate bigWig file from bed file
*/

process BED_TO_BIGWIG {
  label 'bigmem'
  
  input: 
  tuple val(original_vcf), path(bed), val(genome), val(source), val(priorities)
  
  afterScript 'rm all.bed'
  
  shell:
  output_dir = params.output_dir
  output_wig = file(original_vcf).getName().replace(".vcf.gz", ".wig")
  output_bw = file(original_vcf).getName().replace(".vcf.gz", ".bw")
  chrom_sizes = "${moduleDir}/../../common_files/chrom_sizes/${genome}.chrom.sizes"
  
  '''
  bed_to_wig !{bed} !{output_wig}
    
  wigToBigWig -clip -keepAllChromosomes -fixedSummaries \
    !{output_wig} \
    !{chrom_sizes} \
    !{output_bw}
  
  mkdir -p !{output_dir}/!{genome}/!{source}/tracks
  mv !{output_bw} !{output_dir}/!{genome}/!{source}/tracks/
  '''
}
