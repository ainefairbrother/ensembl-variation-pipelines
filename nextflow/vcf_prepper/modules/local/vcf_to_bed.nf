#!/usr/bin/env nextflow

process VCF_TO_BED {
  input: 
  path rank_file
  tuple val(meta), path(vcf) 
  
  output:
  tuple val(meta), path(output_file)
  
  shell:
  output_file = vcf.getName().replace(".vcf.gz", ".bed")
  
  '''
  vcf_to_bed !{vcf} !{output_file} !{rank_file}
    
  rm !{vcf}
  '''
}
