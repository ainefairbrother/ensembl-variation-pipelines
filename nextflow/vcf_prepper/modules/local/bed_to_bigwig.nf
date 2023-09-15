#!/usr/bin/env nextflow

process BED_TO_BIGWIG {
  label 'bigmem'
  
  input: 
  tuple val(meta), path(bed)
  
  output:
  path "variant-${source}-summary.bw"
  
  afterScript 'rm all.bed'
  
  shell:
  source = meta.source
  output_wig = "variant-${source}-summary.wig"
  output_bw = "${meta.genome_tracks_outdir}/variant-${source}-summary.bw"
  chrom_sizes = meta.chrom_sizes
  
  '''
  bed_to_wig !{bed} !{output_wig}
    
  wigToBigWig -clip -keepAllChromosomes -fixedSummaries \
    !{output_wig} \
    !{chrom_sizes} \
    !{output_bw}
    
  ln -sf !{output_bw} "variant-!{source}-summary.bw"
  '''
}
