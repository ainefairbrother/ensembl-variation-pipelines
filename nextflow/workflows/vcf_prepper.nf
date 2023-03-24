/* Workflow to run Ensembl Variation vcf prepper Pipeline
* The Goal of this workflow is to process and annotate VCF files and generate related track files for genome browser
*/

nextflow.enable.dsl=2

import groovy.json.JsonSlurper
import java.io.File

// params default
params.input_config = "${projectDir}/../nf_config/input_sources.json"
params.output_dir = "/nfs/production/flicek/ensembl/variation/new_website"
params.state = "pre"

// params for nextflow-vep
params.vep_config = "${projectDir}/../nf_config/vep.ini"
params.outdir = params.output_dir

// module imports
include { mergeVCF } from "${projectDir}/../nf_modules/merge_vcf.nf"
include { runVEP } from "${projectDir}/../nf_modules/run_vep.nf"
include { generateTracks } from "${projectDir}/../nf_modules/generate_tracks.nf"
include { filterChr } from "${projectDir}/../nf_modules/filter_chr.nf"
include { renameChr } from "${projectDir}/../nf_modules/rename_chr.nf"
include { removeDupIDs } from "${projectDir}/../nf_modules/remove_dup_ids.nf"


log.info 'Starting workflow.....'

def slurper = new JsonSlurper()
params.config = slurper.parse(new File(params.input_config))

workflow {
  genomes = params.config.keySet()
  for (genome in genomes) {
    vcf_files = []
    for (source in params.config.get(genome)){
      vcf_files.add([source.source_name, source.file_location])
    }
    
    // create a directory for this genome in the output directory
    genome_outdir = params.output_dir + "/" + genome
    file(genome_outdir).mkdir()
    
    ch = Channel.fromList(vcf_files)
    
    state = params.state
    //if ( state.equals("pre") ) {
    //  mergeVCF(ch, source_vcf_outdir)
    //  state = "vep"
    //}
    //if ( state.equals("vep") ) {
    //  runVEP(ch, genome_outdir)
    //  state = "post"
    //}
    if( state.equals("post")  ) {
      renameChr(ch, genome_outdir)
      filterChr(renameChr.out.vcfFile)
      removeDupIDs(filterChr.out.vcfFile)
      state = "focus"
    }
    // generateTracks(ch)
  }
}
