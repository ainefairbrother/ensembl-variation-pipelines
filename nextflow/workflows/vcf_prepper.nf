/* Workflow to run Ensembl Variation vcf prepper Pipeline
* The Goal of this workflow is to process and annotate VCF files and generate related track files for genome browser
*/

nextflow.enable.dsl=2

import groovy.json.JsonSlurper
import java.io.File

// params default
params.input_config = "${projectDir}/../nf_config/input_sources.json"
params.output_dir = "/nfs/production/flicek/ensembl/variation/new_website"

// params for nextflow-vep
params.singularity_dir = "/hps/nobackup/flicek/ensembl/variation/snhossain/website/singularity-images"
params.bin_size = 250000
params.remove_patch = 1

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
// tracks
include { readChrVCF } from "${projectDir}/../nf_modules/read_chr_vcf.nf"
include { splitChrVCF } from "${projectDir}/../nf_modules/split_chr_vcf.nf"
include { vcfToBed } from "${projectDir}/../nf_modules/vcf_to_bed.nf"
include { concatBed } from "${projectDir}/../nf_modules/concat_bed.nf"
include { bedToBigBed } from "${projectDir}/../nf_modules/bed_to_bigbed.nf"
include { bedToBigWig } from "${projectDir}/../nf_modules/bed_to_bigwig.nf"
// focus track
include { createFocusTrack } from "${projectDir}/../nf_modules/create_focus_track.nf"


log.info 'Starting workflow.....'

def slurper = new JsonSlurper()
params.config = slurper.parse(new File(params.input_config))

workflow {
  ch_params = [:]
  ch_params['vcf_files'] = []
  ch_params['outdirs'] = []
  ch_params['vep_ini_files'] = []
  
  priorities = [:]

  genomes = params.config.keySet()
  for (genome in genomes) {

    // create a directory for this genome in the output directory
    genome_outdir = params.output_dir + "/" + genome
    file(genome_outdir).mkdir()
    
    priorities[genome] = [:]
  
    for (source in params.config.get(genome)){
      source_name = source.source_name
    
      // create a directory for this genome in the output directory
      vep_outdir = genome_outdir + "/" + source_name.replace(" ", "_") + "/vcfs"
      file(vep_outdir).mkdir()
      
      // check if index file exists
      index_file = ""
      if (!file(source.file_location + ".tbi").exists() && !file(source.file_location + ".csi").exists()){
        exit 1, "index file does not exist for - " + source.file_location
      }
      
      // vep config file
      vep_ini = "${projectDir}/../nf_config/vep_ini/${genome}.ini"
      
      ch_params['vcf_files'].add(source.file_location)
      ch_params['outdirs'].add(vep_outdir)
      ch_params['vep_ini_files'].add(vep_ini)
      
      // priority of a source in a genome - used in focus track creation
      priorities[genome][source_name] = source.priority ? source.priority : 100
    }
  }

  vcf_files = Channel.fromList(ch_params['vcf_files'])
  outdirs = Channel.fromList(ch_params['outdirs'])
  vep_configs = Channel.fromList(ch_params['vep_ini_files'])

  //mergeVCF(ch, source_vcf_outdir)
  
  vep(vcf_files, vep_configs, outdirs)
  
  renameChr(
      vep.out,
      vep.out.map{ file(it).getParent().getParent().getParent().getSimpleName() },
      vep.out.map{ file(it).getParent().getParent().getSimpleName() },
      priorities
    )
  removeDupIDs(renameChr.out)
  indexVCF(removeDupIDs.out)
  
  readChrVCF(indexVCF.out)
  splitChrVCF(readChrVCF.out.transpose())
  vcfToBed(splitChrVCF.out.transpose())
  concatBed(vcfToBed.out.groupTuple(by: [0, 2, 3]))
    
  bedToBigBed(concatBed.out)
  bedToBigWig(concatBed.out)
  
  createFocusTrack(concatBed.out.groupTuple(by: [2]))  
}
