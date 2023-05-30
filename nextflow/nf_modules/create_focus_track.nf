#!/usr/bin/env nextflow

/*
* This script create combined focus track files from multiple sources
*/

process createFocusTrack {
  label 'bigmem'

  input:
  tuple val(original_vcfs), path(bed_files), val(genome), val(sources), val(priorities)
  
  afterScript 'rm all.bed'
  
  shell:
  output_dir = params.output_dir
  output_bed = genome + ".bed"
  output_bb = genome + ".bb"
  output_wig = genome + ".wig"
  output_bw = genome + ".bw"
  
  '''
  priorities="!{priorities}"
  priorities=(${priorities//[,\\[\\]]/})
  bed_files=(!{bed_files})
  
  # we need to order the bed files with priority
  total_idx=${#priorities[@]}
  let "total_idx--"
  for idx in $(seq 1 ${total_idx});
  do
    prev_idx=$((idx-1))
    if [[ ${priorities[${idx}]} < ${priorities[${prev_idx}]} ]]; then
      temp=${priorities[${idx}]}
      priorities[${idx}]=${priorities[${prev_idx}]}
      priorities[${prev_idx}]=${temp}
      
      temp=${bed_files[${idx}]}
      bed_files[${idx}]=${bed_files[${prev_idx}]}
      bed_files[${prev_idx}]=${temp}
    fi
  done
  
  !{projectDir}/../../bin/merge_bed \
    all.bed \
    ${bed_files[@]} \
    
  LC_COLLATE=C sort -S1G -k1,1 -k2,2n all.bed > !{output_bed}
    
  chrom_sizes=!{projectDir}/../nf_config/chrom_sizes/!{genome}.chrom.sizes
  
  bedToBigBed -type=bed3+6 !{output_bed} ${chrom_sizes} !{output_bb}
  
  !{projectDir}/../../bin/bed_to_wig \
    !{output_bed} \
    !{output_wig}
    
  wigToBigWig -clip -keepAllChromosomes -fixedSummaries \
    !{output_wig} \
    ${chrom_sizes} \
    !{output_bw}
    
  mkdir -p !{output_dir}/!{genome}/focus/tracks
  mv !{output_bb} !{output_dir}/!{genome}/focus/tracks/
  mv !{output_bw} !{output_dir}/!{genome}/focus/tracks/
  '''
}