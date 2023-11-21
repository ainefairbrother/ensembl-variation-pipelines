//
// This file holds several functions specific to the main.nf workflow in the vcf_prepper pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // TODO: Validate parameters and print summary to screen
    // Create output directory structure
    //
    public static void initialise(workflow, params, log) {
        
        // Check proper directory paths have been provided has been provided
        if (!params.output_dir) {
            log.error "Please provide a directory path to create output files e.g. '--output_dir OUTPUT_DIR'"
            System.exit(1)
        }
        
        if (!params.temp_dir) {
            log.error "Please provide a directory path to create tmp files e.g. '--temp_dir TMP_DIR'"
            System.exit(1)
        }
        
        if( !Nextflow.file(params.output_dir).exists() ) {
            log.warn "output directory path does not exist - ${params.output_dir}, creating ..."
            Nextflow.file(params.output_dir).mkdirs()
        }
        
        if( !Nextflow.file(params.temp_dir).exists() ) {
            log.warn "tmp directory path does not exist - ${params.temp_dir}, creating ..."
            Nextflow.file(params.temp_dir).mkdirs()
        }
        
        // Create output directory structure
        Nextflow.file("${params.output_dir}/api").mkdirs()
        Nextflow.file("${params.output_dir}/tracks").mkdirs()
    }
}