//
// Workflow to run Ensembl Variation vcf prepper Pipeline
// The Goal of this workflow is to process and annotate VCF files and generate related track files for genome browser
//

import groovy.json.JsonSlurper
import java.io.File

def slurper = new JsonSlurper()
params.config = slurper.parse(new File(params.input_config))

include { CREATE_RANK_FILE } from "../modules/local/create_rank_file.nf"
include { PREPARE_GENOME } from "../subworkflows/local/prepare_genome.nf"
include { UPDATE_FIELDS } from "../modules/local/update_fields.nf"
include { REMOVE_VARIANTS } from "../modules/local/remove_variants.nf"
include { RUN_VEP } from "../subworkflows/local/run_vep.nf"
include { SPLIT_VCF } from "../subworkflows/local/split_vcf.nf"
include { VCF_TO_BED } from "../modules/local/vcf_to_bed.nf"
include { CONCAT_BEDS } from "../modules/local/concat_beds.nf"
include { BED_TO_BIGBED } from "../modules/local/bed_to_bigbed.nf"
include { BED_TO_BIGWIG } from "../modules/local/bed_to_bigwig.nf"

def parse_config (config) {
  input_set = []
  
  genomes = config.keySet()
  for (genome in genomes) {
    for (source_data in params.config.get(genome)) {
      vcf = source_data.file_location
      
      meta = [:]
      meta.genome = genome
      meta.genome_uuid = source_data.genome_uuid
      meta.species = source_data.species
      meta.assembly = source_data.assembly
      meta.source = source_data.source_name.replace(" ", "_")
      meta.file_type = source_data.file_type
      
      input_set.add([meta, vcf])
    }  
  }
  
  return input_set
}

workflow VCF_PREPPER {
  if (params.skip_vep && params.skip_tracks) {
    exit 0, "Skipping VEP and track file generation, nothing to do ..."
  }
  
  input_set = parse_config(params.config)
  ch_input = Channel.fromList( input_set )
  
  // setup
  CREATE_RANK_FILE( params.rank_file )
  PREPARE_GENOME( ch_input )
    
  // api files
  if (!params.skip_vep) {
    // pre-process
    UPDATE_FIELDS( PREPARE_GENOME.out )
    REMOVE_VARIANTS( UPDATE_FIELDS.out )
    
    // run vep
    vep = RUN_VEP( REMOVE_VARIANTS.out )
    vep.view()
    // post-process
    vep
    .map {
      meta, vcf, vcf_index ->
        // TODO: when we have multiple source per genome we need to delete source specific files
        new_vcf = "${meta.genome_api_outdir}/variation.vcf.gz"
        new_vcf_index = "${meta.genome_api_outdir}/variation.vcf.gz.${meta.index_type}"
        
        // in -resume vcf and vcf_index may not exists as already renamed
        // moveTo instead of renameTo - in -resume dest file may exists from previous run
        if ( file(vcf).exists() && file(vcf_index).exists() ) {
          file(vcf).moveTo(new_vcf)
          file(vcf_index).moveTo(new_vcf_index)
        }
        [meta, new_vcf, new_vcf_index]
    }.set { ch_tracks }
  }
  else {
    ch_tracks = PREPARE_GENOME.out
  }
  
  // track files
  if (!params.skip_tracks) {
    // create bed from VCF
    // TODO: vcf_to_bed maybe faster without SPLIT_VCF - needs benchmarking
    SPLIT_VCF( ch_tracks )
    VCF_TO_BED( CREATE_RANK_FILE.out, SPLIT_VCF.out.transpose() )
    CONCAT_BEDS( VCF_TO_BED.out.groupTuple() )

    // create source tracks
    // TODO: remove symlink creation for focus track when we have multiple source
    BED_TO_BIGBED( CONCAT_BEDS.out )
    BED_TO_BIGWIG( CONCAT_BEDS.out )
  }
}
