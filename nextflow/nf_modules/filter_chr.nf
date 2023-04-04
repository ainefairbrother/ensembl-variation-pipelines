#!/usr/bin/env nextflow

/*
* This script filter out variants from some specific seq region (patches, contigs, etc.)
*/

process filterChr {
  input:
  tuple val(original), path(vcfFile), path(indexFile)
  
  output:
  tuple val(original), path("filtered-${original}-*.vcf"), emit: file
  
  shell:
  '''
  vcf_file=!{vcfFile}
  output_file=filtered-!{original}-${vcf_file/%.gz/}
  
  # seq regions that we want to keep
  keep_regions=$(tabix !{vcfFile} -l | grep -v -E '(CTG|PATCH|TEST)' |  paste -sd,)
  
  # only keep the regions that we want to keep
  bcftools view -r ${keep_regions} !{vcfFile} -Ov -o ${output_file}
  '''
}
