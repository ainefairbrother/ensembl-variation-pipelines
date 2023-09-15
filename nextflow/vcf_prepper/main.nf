include { VCF_PREPPER } from "./workflows/vcf_prepper.nf"

WorkflowMain.initialise(workflow, params, log)

workflow {
  VCF_PREPPER()
}