/* Workflow to run Ensembl Variation vcf prepper Pipeline
* The Goal of this workflow is to process and annotate VCF files and generate related track files for genome browser
*/

nextflow.enable.dsl=2

import groovy.json.JsonSlurper

// params default

params.input_config = "${projectDir}/../nf_config/input_sources.json"

// module imports
include { generateTracks } from "${projectDir}/../nf_modules/generate_tracks.nf"


log.info 'Starting workflow.....'

def slurper = new JsonSlurper()
params.config = slurper.parse(new File(params.input_config))

workflow {
  genomes = params.config.keySet()
  vcf_files = []
  for (genome in genomes){
    for (source in params.config.get(genome)){
      vcf_files.add([source.source_name, source.file_location])
    }
  }
  
  ch = Channel.fromList(vcf_files)
  generateTracks(ch)
}