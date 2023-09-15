//
// run VEP annotation
//

repo_dir = params.repo_dir

include { INDEX_VCF } from "../../modules/local/index_vcf.nf"
include { vep } from "${repo_dir}/ensembl-vep/nextflow/workflows/run_vep.nf"

workflow RUN_VEP {
  take:
    input
  
  main:
  // create index file at the exact location where input vcf file is as nextflow-vep requires as such
  input
  .map {
    meta, vcf ->
      vcf_fullpath = vcf.toString()
      [meta, vcf_fullpath]
  }
  .set { ch_index_vcf }
  INDEX_VCF( ch_index_vcf )
  
  INDEX_VCF.out
  .multiMap {
    meta, vcf, vcf_index ->
      inputs:
        // set vcf file name as <genome>-<source> so that we can track meta later on
        new_vcf = file(vcf).getParent() + "/${meta.genome}-${meta.source}.vcf.gz"
        new_vcf_index = "${new_vcf}.${meta.index_type}"
        
        // moveTo instead of renameTo because in -resume dest file may already exists in work dir
        file(vcf).moveTo(new_vcf)
        file(vcf_index).moveTo(new_vcf_index)
        
        // TODO: for -resume abitlity create symlink
        //file(new_vcf).mklink(vcf.toString())
        //file(new_vcf_index).mklink(vcf_index.toString())
        
        [new_vcf, meta.vep_config] 
      output_dir: meta.genome_api_outdir
      one_to_many: 0
  }
  .set { ch_vep }
  vep( ch_vep.inputs, ch_vep.output_dir, ch_vep.one_to_many )
  
  input
  .map {
    meta, vcf ->
      // tag here is the output vcf file from nextflow-vep
      tag = "${meta.genome_api_outdir}/${meta.genome}-${meta.source}_VEP.vcf.gz"    
      [tag, meta]
  }
  .join ( vep.out, failOnDuplicate: true )
  .map {
    tag, meta ->
      vcf = tag
      vcf_index = "${tag}.${meta.index_type}"
      
      [meta, vcf, vcf_index]
  }.set { ch_post_vep }
  
  emit:
    ch_post_vep
}