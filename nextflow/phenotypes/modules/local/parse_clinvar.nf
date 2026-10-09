process PARSE_CLINVAR {
    tag "${rcv_xml.name} (${assembly})"

    input:
    path rcv_xml
    path vcv_xml
    path python_source
    val assembly

    output:
    path 'clinvar_output.jsonl', emit: records
    path 'clinvar_summary.json', emit: summary
    path 'warnings.txt', emit: warnings
    path 'rejections.txt', emit: rejections

    script:
    """
    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="${python_source}" \
    "${projectDir}/.venv/bin/python" -m phenotypes.clinvar.rcv_parser \
        --rcv-input "${rcv_xml}" \
        --vcv-input "${vcv_xml}" \
        --assembly "${assembly}" \
        --records clinvar_output.jsonl \
        --summary clinvar_summary.json
    """
}
