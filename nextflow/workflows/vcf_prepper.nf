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

// module imports
include { mergeVCF } from "${projectDir}/../nf_modules/merge_vcf.nf"
include { generateTracks } from "${projectDir}/../nf_modules/generate_tracks.nf"


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
    outdir = file(genome_outdir)
    outdir.mkdir()
    
    ch = Channel.fromList(vcf_files)
    
    if ( params.state.equals("pre") ){
      mergeVCF(ch, genome_outdir)
      params.state = "vep"
    }
    // runVEP(ch)
    // generateTracks(ch)
  }
}