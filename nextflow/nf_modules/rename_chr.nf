#!/usr/bin/env nextflow

/*
* This script rename the seq region synonyms in the VCF file
*/

process renameChr {
  input:
  val input_file
  val genome
  val source
  val priorities
  
  output:
  tuple val(output_dir), val(prefix), val(genome), val(source), val(priority), env(index_type)
  
  shell:
  priority = priorities[genome][source]
  output_dir = params.output_dir + "${genome}/${source}/vcfs/"
  prefix = file(input_file).getSimpleName()
  
  '''
  output_file=!{output_dir}/!{prefix}_renamed_VEP.vcf.gz
  synonym_file=!{projectDir}/../nf_config/synonyms/!{genome}.txt
  
  # rename chr synonyms
  bcftools annotate --no-version --force --rename-chrs ${synonym_file} !{input_file} -Oz -o ${output_file}
  
  # for next job create parameters
  if [[ -f !{input_file}.tbi ]]; then
    index_type=tbi
  else
    index_type=csi
  fi
  '''
}
