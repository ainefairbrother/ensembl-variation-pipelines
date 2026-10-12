process PARSE_CLINVAR {
    tag "ClinVar ${processing_context.species_production_name} ${processing_context.assembly}"

    // Resolve the path per task, after the processing context has been supplied.
    publishDir path: { "${params.outdir}/clinvar/${processing_context.species_production_name}/${processing_context.assembly}" },
        mode: 'copy', overwrite: true

    input:
    path rcv_xml
    path vcv_xml
    path python_source
    val python_source_contents
    val processing_context

    output:
    tuple val(processing_context), path('clinvar_output.jsonl'), emit: records
    tuple val(processing_context), path('clinvar_summary.json'), emit: summary
    tuple val(processing_context), path('warnings.txt'), emit: warnings
    tuple val(processing_context), path('rejections.txt'), emit: rejections

    script:
    // Configured fields already match GenomeContext. Pass them unchanged,
    // quoted as one literal shell argument, including any apostrophes.
    def context_argument = "'" + groovy.json.JsonOutput.toJson(processing_context).replace("'", "'\"'\"'") + "'"
    """
    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="${python_source}" \
    "${projectDir}/.venv/bin/python" -m phenotypes.clinvar.rcv_parser \
        --rcv-input "${rcv_xml}" \
        --vcv-input "${vcv_xml}" \
        --genome-context ${context_argument} \
        --pipeline-run-date "${workflow.start.format(java.time.format.DateTimeFormatter.ISO_LOCAL_DATE)}" \
        --pipeline-directory "${projectDir}" \
        --records clinvar_output.jsonl \
        --summary clinvar_summary.json
    """
}
