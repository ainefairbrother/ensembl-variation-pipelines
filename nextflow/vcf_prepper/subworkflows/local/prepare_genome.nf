//
// Prepare genome related files
//

include { GENERATE_CHROM_SIZES } from "../../modules/local/generate_chrom_sizes.nf"
include { GENERATE_VEP_CONFIG } from "../../modules/local/generate_vep_config.nf"
include { GENERATE_SYNONYM_FILE } from "../../modules/local/generate_synonym_file.nf"
include { PROCESS_CACHE } from "../../modules/local/process_cache.nf"
include { PROCESS_FASTA } from "../../modules/local/process_fasta.nf"
include { PROCESS_CONSERVATION_DATA } from "../../modules/local/process_conservation_data.nf"
include { DOWNLOAD_SOURCE } from "../../modules/local/download_source.nf"

workflow PREPARE_GENOME {
  take:
    input
  
  main:
    input
    .map {
      meta, vcf ->
        genome_temp_dir = "${params.temp_dir}/${meta.genome_uuid}"
        file(genome_temp_dir).mkdirs()
        
        synonym_file = "${genome_temp_dir}/${meta.genome}.synonyms"
        vep_config = "${genome_temp_dir}/${meta.genome}.ini"
        chrom_sizes = "${genome_temp_dir}/${meta.genome}.chrom.sizes"
        
        genome_api_outdir = "${params.output_dir}/api/${meta.genome_uuid}"
        file(genome_api_outdir).mkdirs()
        genome_tracks_outdir = "${params.output_dir}/tracks/${meta.genome_uuid}"
        file(genome_tracks_outdir).mkdirs()
        
        cache_dir = if meta.species == "homo_sapiens" ? params.cache_dir : genome_temp_dir
        fasta_dir = if meta.species == "homo_sapiens" ? params.fasta_dir : genome_temp_dir
        conservation_data_dir = if meta.species == "homo_sapiens" ? params.conservation_data_dir : genome_temp_dir
        
        [ meta + [
            synonym_file: synonym_file,
            vep_config: vep_config,
            chrom_sizes: chrom_sizes,
            genome_temp_dir: genome_temp_dir,
            genome_api_outdir: genome_api_outdir,
            genome_tracks_outdir: genome_tracks_outdir,
            cache_dir: cache_dir,
            fasta_dir: fasta_dir,
            conservation_data_dir: conservation_data_dir
          ], vcf
        ]
    }.set { ch_prepare_genome }
    
    // TODO: run this only once per genome - currently they run for all sources in a genome
    ch_vep_config_done = GENERATE_VEP_CONFIG( ch_prepare_genome.map { meta, vcf -> meta } )
    ch_synonym_file_done = GENERATE_SYNONYM_FILE( ch_prepare_genome.map { meta, vcf -> meta } )
    ch_chrom_sizes_done = GENERATE_CHROM_SIZES( ch_prepare_genome.map { meta, vcf -> meta } )
    ch_processed_cache = PROCESS_CACHE( ch_prepare_genome.map { meta, vcf -> meta } )
    ch_processed_fasta = PROCESS_FASTA( ch_prepare_genome.map { meta, vcf -> meta } )
    ch_processed_conservation = PROCESS_CONSERVATION_DATA( ch_prepare_genome.map { meta, vcf -> meta } )
    
    // we join channels to only create DAG edges
    ch_prepare_genome
    .map {
      meta, vcf ->
        // tag here is the genome
        tag = meta.genome 
          
        [tag, meta, vcf]
    }
    .join ( ch_vep_config_done )
    .join ( ch_synonym_file_done )
    .join ( ch_chrom_sizes_done )
    .join ( ch_processed_cache )
    .join ( ch_processed_fasta )
    .join ( ch_processed_conservation )
    .map {
      tag, meta, vcf ->
        [meta, vcf]
    }
    .set { ch_prepare_source }
    
    DOWNLOAD_SOURCE( ch_prepare_input )
    
    emit:
      DOWNLOAD_SOURCE.out
}