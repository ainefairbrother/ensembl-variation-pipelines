#!/usr/bin/env nextflow

process BED_TO_BIGBED {
  input: 
  tuple val(meta), path(bed)
  
  output:
  path "variant-${source}-details.bb"
  
  shell:
  source = meta.source
  output_bb = "${meta.genome_tracks_outdir}/variant-${source}-details.bb"
  chrom_sizes = meta.chrom_sizes
  
  '''
  bedToBigBed -type=bed3+6 !{bed} !{chrom_sizes} !{output_bb}
  ln -sf !{output_bb} "variant-!{source}-details.bb"
  
  # temp: for one source we create symlink for focus
  cd !{meta.genome_tracks_outdir}
  ln -sf variant-!{source}-details.bb variant-details.bb
  '''
}
