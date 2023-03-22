#!/usr/bin/env nextflow

/*
* This script run nextflow-vep on VCF file
*/

process runVEP {
  input: 
  tuple val(source), val(vcfFile)
  val genome_dir
  
  shell:
  '''
  vcf_file=!{vcfFile}
  
  # remove any whitespace in source name
  source="!{source}"
  source=${source// /_}
  
  # create a source dir under output dir
  source_dir=!{genome_dir}/${source}
  mkdir -p $source_dir
  
  # create output dir
  output_dir=${source_dir}/vcfs
  mkdir -p $output_dir
  
  # format output file prefix
  output_file_prefix=$(basename ${vcf_file/.vcf.gz/_VEP})
  output_file_prefix=${output_file_prefix/_no_VEP/}
  
  repositories="/hps/software/users/ensembl/repositories/${USER}"
  
  module load openjdk-1.8.0_265-b01-gcc-9.3.0-w3gaafy
  
  cd ${output_dir}
  nextflow -C ${repositories}/ensembl-vep/nextflow/nf_config/nextflow.config \\
    run ${repositories}/ensembl-vep/nextflow/workflows/run_vep.nf \\
    --vcf ${vcf_file} \\
    --vep_config !{projectDir}/../nf_config/vep.ini \\
    --outdir ${output_dir} \\
    --output_prefix ${output_file_prefix} \\
    --singularity_dir /hps/nobackup/flicek/ensembl/variation/snhossain/website/singularity-images
    -profile lsf \\
    --bin_size 1000000 \\
    -resume \\
    -w work_${source}
  '''
}