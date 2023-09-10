/* Workflow to run Ensembl Variation vcf prepper Pipeline
* The Goal of this workflow is to process and annotate VCF files and generate related track files for genome browser
*/

nextflow.enable.dsl=2

import groovy.json.JsonSlurper
import java.io.File

// params default
params.input_config = "${projectDir}/common_files/input_sources.json"
params.output_dir = "/nfs/production/flicek/ensembl/variation/new_website"
params.singularity_dir = "/hps/nobackup/flicek/ensembl/variation/snhossain/website/singularity-images"
params.bin_size = 250000
params.remove_nonunique_ids = 1
params.remove_patch_regions = 1
params.skip_vep = 0
params.skip_tracks = 0
params.skip_create_config = 0
params.rename_clinvar_ids = 1
params.ini_file = "${projectDir}/common_files/DEFAULT.ini"
params.rank_file = "${projectDir}/common_files/variation_consequnce_rank.json"
params.version = 108

// module imports
repo_dir = "/hps/software/users/ensembl/repositories/${USER}"
// pre
include { CREATE_RANK_FILE } from "./modules/local/create_rank_file.nf"
include { CREATE_CONFIGS } from "./modules/local/create_configs.nf"
include { MERGE_VCFS } from "./modules/local/merge_vcfs.nf"
// vep
include { vep } from "${repo_dir}/ensembl-vep/nextflow/workflows/run_vep.nf"
// post
include { RENAME_CHR } from "./modules/local/rename_chr.nf"
include { REMOVE_VARIANTS } from "./modules/local/remove_variants.nf"
include { RENAME_CLINVAR_IDS } from "./modules/local/rename_clinvar_ids.nf"
include { INDEX_VCF } from "./modules/local/index_vcf.nf"
// tracks
include { READ_CHR } from "./modules/local/read_chr.nf"
include { SPLIT_VCF } from "./modules/local/split_vcf.nf"
include { VCF_TO_BED } from "./modules/local/vcf_to_bed.nf"
include { CONCAT_BEDS } from "./modules/local/concat_beds.nf"
include { BED_TO_BIGBED } from "./modules/local/bed_to_bigbed.nf"
include { BED_TO_BIGWIG } from "./modules/local/bed_to_bigwig.nf"
// focus track
include { CREATE_FOCUS_TRACK } from "./modules/local/create_focus_track.nf"


log.info 'Starting workflow.....'

def slurper = new JsonSlurper()
params.config = slurper.parse(new File(params.input_config))

workflow {
  if (params.skip_vep && params.skip_tracks){
    exit 0, "Skipping VEP and track file generation, nothing to do ..."
  }

  ch_params = [:]
  ch_params['vcf_files'] = []
  ch_params['outdirs'] = []
  ch_params['vep_ini_files'] = []
  ch_params['genomes'] = []
  ch_params['sources'] = []
  
  priorities = [:]

  genomes = params.config.keySet()
  for (genome in genomes) {

    // create a directory for this genome in the output directory
    genome_outdir = params.output_dir + "/" + genome
    file(genome_outdir).mkdir()
    
    priorities[genome] = [:]
  
    for (source in params.config.get(genome)){
      source_name = source.source_name
      source_name = source_name.replace(" ", "_") 
    
      // create a directory for this genome in the output directory
      vep_outdir = genome_outdir + "/" + source_name + "/vcfs"
      file(vep_outdir).mkdirs()
      
      // check if index file exists
      index_file = ""
      if (!file(source.file_location + ".tbi").exists() && !file(source.file_location + ".csi").exists()){
        exit 1, "index file does not exist for - " + source.file_location
      }
      
      // vep config file
      vep_ini = "${projectDir}/common_files/vep_ini/${genome}.ini"
      
      ch_params['vcf_files'].add(source.file_location)
      ch_params['outdirs'].add(vep_outdir)
      ch_params['vep_ini_files'].add(vep_ini)
      ch_params['genomes'].add(genome)
      ch_params['sources'].add(source_name)
      
      // priority of a source in a genome - used in focus track creation
      priorities[genome][source_name] = source.priority ? source.priority : 100
    }
  }

  // set up Channels
  vcf_files = Channel.fromList(ch_params['vcf_files'])
  outdirs = Channel.fromList(ch_params['outdirs'])
  vep_configs = Channel.fromList(ch_params['vep_ini_files'])
  genomes = Channel.fromList(ch_params['genomes'])
  sources = Channel.fromList(ch_params['sources'])

  // create config files
  CREATE_RANK_FILE(params.rank_file)
  CREATE_CONFIGS(CREATE_RANK_FILE.out, genomes)
  
  //MERGE_VCFS(ch, source_vcf_outdir)
  
  // run VEP on input VCF file
  if (!params.skip_vep) {
    vep(vcf_files, vep_configs, outdirs)
  }
  
  // create track files from VEPed VCF file
  // first, post process the VEPed file
  if (!params.skip_tracks){
    if(!params.skip_vep){
      RENAME_CHR(
        vep.out,
        vep.out.map{ file(it).getParent().getParent().getParent().getSimpleName() },
        vep.out.map{ file(it).getParent().getParent().getSimpleName() },
        priorities
      )
    }
    else {
      RENAME_CHR(
        vcf_files,
        genomes,
        sources,
        priorities
      )
    }
    REMOVE_VARIANTS(RENAME_CHR.out)
    RENAME_CLINVAR_IDS(REMOVE_VARIANTS.out)
    INDEX_VCF(RENAME_CLINVAR_IDS.out)
    
    // then, we create bed files from the each VEPed VCF file
    READ_CHR(INDEX_VCF.out)
    SPLIT_VCF(READ_CHR.out.transpose())
    VCF_TO_BED(CREATE_CONFIGS.out.collect(), SPLIT_VCF.out.transpose())
    CONCAT_BEDS(VCF_TO_BED.out.groupTuple(by: [0, 2, 3, 4]))

    // then, we create bigBed and bigWig files from each bed file
    BED_TO_BIGBED(CONCAT_BEDS.out)
    BED_TO_BIGWIG(CONCAT_BEDS.out)

    // in the same time, we also create the focus track files from concatenated bed file
    CREATE_FOCUS_TRACK(CONCAT_BEDS.out.groupTuple(by: [2]))
  } 
}
