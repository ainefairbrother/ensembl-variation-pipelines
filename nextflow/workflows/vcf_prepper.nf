/* Workflow to run Ensembl Variation vcf prepper Pipeline
* The Goal of this workflow is to process and annotate VCF files and generate related track files for genome browser
*/

nextflow.enable.dsl=2

import groovy.json.JsonSlurper
import java.io.File

// params default
params.input_config = "${projectDir}/../nf_config/input_sources.json"
//params.output_dir = "/nfs/production/flicek/ensembl/variation/new_website"
params.output_dir = "/hps/nobackup/flicek/ensembl/variation/snhossain/website/test_v108"
params.state = "pre"

// params for nextflow-vep
params.singularity_dir = "/hps/nobackup/flicek/ensembl/variation/snhossain/website/singularity-images"
params.bin_size = 100

// module imports
repo_dir = "/hps/software/users/ensembl/repositories/${USER}"
// pre
include { mergeVCF } from "${projectDir}/../nf_modules/merge_vcf.nf"
// vep
include { vep } from "${repo_dir}/ensembl-vep/nextflow/workflows/run_vep.nf"
// post
include { renameChr } from "${projectDir}/../nf_modules/rename_chr.nf"
include { removeDupIDs } from "${projectDir}/../nf_modules/remove_dup_ids.nf"
include { indexVCF } from "${projectDir}/../nf_modules/index_vcf.nf"


//include { checkVCF; checkVCF } from "${repo_dir}/ensembl-vep/nextflow/nf_modules/check_VCF.nf"
//include { readVCF } from "${repo_dir}/ensembl-vep/nextflow/nf_modules/read_VCF.nf"
//include { splitVCF; splitVCF } from "${repo_dir}/ensembl-vep/nextflow/nf_modules/split_VCF.nf"
//include { mergeVCF } from "${repo_dir}/ensembl-vep/nextflow/nf_modules/merge_VCF.nf"

// tracks
include { readChrVCF } from "${projectDir}/../nf_modules/read_chr_vcf.nf"
include { splitChrVCF } from "${projectDir}/../nf_modules/split_chr_vcf.nf"
include { vcfToBed } from "${projectDir}/../nf_modules/vcf_to_bed.nf"
include { bedToBigBed } from "${projectDir}/../nf_modules/bed_to_bigbed.nf"
include { createFocusVCF } from "${projectDir}/../nf_modules/create_focus_vcf.nf"


log.info 'Starting workflow.....'

def slurper = new JsonSlurper()
params.config = slurper.parse(new File(params.input_config))

workflow {
  ch_params = [:]
  ch_params['vcf_files'] = []
  ch_params['index_files'] = []
  ch_params['outdirs'] = []
  ch_params['input_prefixes'] = []
  ch_params['vep_ini_files'] = []
  ch_params['synonyms'] = []

  genomes = params.config.keySet()
  for (genome in genomes) {

    // create a directory for this genome in the output directory
    genome_outdir = params.output_dir + "/" + genome
    file(genome_outdir).mkdir()
  
    for (source in params.config.get(genome)){
      source_name = source.source_name
    
      // create a directory for this genome in the output directory
      vep_outdir = genome_outdir + "/" + source_name.replace(" ", "_") + "/vcfs"
      file(vep_outdir).mkdir()
      
      // check if index file exists
      index_file = ""
      if (file(source.file_location + ".tbi").exists()){
        index_file = file(source.file_location + ".tbi")
      }
      else if (file(source.file_location + ".csi").exists()){
        index_file = file(source.file_location + ".csi")
      }
      else{
        exit 1, "index file does not exist for - " + source.file_location
      }
      
      // output file prefix for VEP output file
      prefix = file(source.file_location).getSimpleName()
      
      // vep config file
      vep_ini = "${projectDir}/../nf_config/vep_ini/${genome}.ini"
      
      // a tsv file containing chrom names and their synonyms (to be used by bcftools to rename synonyms)
      synonyms = "${projectDir}/../nf_config/synonyms/${genome}.txt"
      
      ch_params['vcf_files'].add(source.file_location)
      ch_params['index_files'].add(index_file)
      ch_params['outdirs'].add(vep_outdir)
      ch_params['input_prefixes'].add(prefix)
      ch_params['vep_ini_files'].add(vep_ini)
      ch_params['synonyms'].add(synonyms)
    }
  }
  
  vcf_files = Channel.fromList(ch_params['vcf_files'])
  index_files = Channel.fromList(ch_params['index_files'])
  outdirs = Channel.fromList(ch_params['outdirs'])
  input_prefixes = Channel.fromList(ch_params['input_prefixes'])
  vep_configs = Channel.fromList(ch_params['vep_ini_files'])
  synonyms = Channel.fromList(ch_params['synonyms'])

  state = params.state

  //if ( state.equals("pre") ) {
  //  mergeVCF(ch, source_vcf_outdir)
  //  state = "vep"
  //}
  if ( state.equals("vep") ) {
    // for multiple vcf this may not be working 
    vep(vcf_files, vep_configs, outdirs)
    
    state = "post"
  }
  if( state.equals("post")  ) {
    renameChr(vep.out, vep.out
      .map{ file(it)
        .getParent()
        .getParent()
        .getParent()
        .getSimpleName() }
    )
    removeDupIDs(renameChr.out)
    indexVCF(removeDupIDs.out)
    
    state = "tracks"
  }
  
  if( state.equals("tracks")  ) {
    readChrVCF(indexVCF.out)
    splitChrVCF(readChrVCF.out.transpose())
    vcfToBed(splitChrVCF.out.transpose())
    
    bedToBigBed(vcfToBed.out.groupTuple())
    //bedToBigWig(vcfToBed.out.groupTuple())
    
    state = "focus"
  }
  if( state.equals("focus")  ) {
    //createFocusVCF(removeDupIDs.out.vcfFile.collect(), removeDupIDs.out.indexFile.collect(), genome_outdir)  
  }
}
