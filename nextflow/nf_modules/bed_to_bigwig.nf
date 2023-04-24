#!/usr/bin/env nextflow

/*
* This script generate bigWig file from bed file
*/

process bedToBigWig {
  input: 
  tuple val(original_vcf), path(bed), val(genome), val(source)
  
  afterScript 'rm all.bed'
  
  shell:
  output_dir = params.output_dir
  output_wig = file(original_vcf).getSimpleName() + ".wig"
  output_bw = file(original_vcf).getSimpleName() + ".bw"
  
  '''
  chrom_sizes=!{projectDir}/../nf_config/chrom_sizes/!{genome}.chrom.sizes
  
  !{projectDir}/../../bin/bed_to_wig \
    !{bed} \
    !{output_wig}
    
  wigToBigWig -clip -keepAllChromosomes -fixedSummaries \
    !{output_wig} \
    ${chrom_sizes} \
    !{output_bw}
  
  mkdir -p !{output_dir}/!{genome}/!{source}/tracks
  mv !{output_bw} !{output_dir}/!{genome}/!{source}/tracks/
  '''
}
