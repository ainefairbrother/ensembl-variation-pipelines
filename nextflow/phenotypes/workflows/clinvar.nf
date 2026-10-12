include { DOWNLOAD_CLINVAR } from '../modules/local/download_clinvar'
include { PARSE_CLINVAR } from '../modules/local/parse_clinvar'

def single_xml_input(pattern, parameter_name) {
    // Collect before selecting: .first() would silently discard extra matches.
    // The resulting value channel can be reused for every configured context.
    return channel.fromPath(pattern, type: 'file', checkIfExists: true)
        .toList()
        .map { paths ->
            if (paths.size() != 1) {
                error("${parameter_name} must match exactly one input file; matched ${paths.size()}.")
            }
            return paths[0]
        }
}

workflow CLINVAR {
    take:
    source_settings

    main:
    releases_to_download = []
    if (!source_settings.inputs.rcv) {
        releases_to_download.add('RCV')
    }
    if (!source_settings.inputs.vcv) {
        releases_to_download.add('VCV')
    }
    if (releases_to_download) {
        DOWNLOAD_CLINVAR(channel.fromList(releases_to_download))
    }

    clinvar_rcv = source_settings.inputs.rcv ? single_xml_input(source_settings.inputs.rcv, 'sources.clinvar.inputs.rcv') :
        DOWNLOAD_CLINVAR.out.xml.filter { release_type, xml -> release_type == 'RCV' }.map { release_type, xml -> xml }.first()
    clinvar_vcv = source_settings.inputs.vcv ? single_xml_input(source_settings.inputs.vcv, 'sources.clinvar.inputs.vcv') :
        DOWNLOAD_CLINVAR.out.xml.filter { release_type, xml -> release_type == 'VCV' }.map { release_type, xml -> xml }.first()

    // Stage the package with its import layout. Include each Python filename
    // and its contents in the cache key; directory metadata misses nested edits.
    python_source = file("${projectDir}/src", checkIfExists: true)
    python_source_contents = channel.fromPath("${projectDir}/src/**/*.py", checkIfExists: true)
        .map { source_file -> tuple(python_source.relativize(source_file).toString(), source_file.text) }
        .toSortedList { left, right -> left[0] <=> right[0] }
    processing_contexts = channel.fromList(source_settings.contexts)
    // Reuse the same source files for every context, rather than consuming
    // one XML channel item for only the first context. Downloads run once.
    PARSE_CLINVAR(clinvar_rcv, clinvar_vcv, python_source, python_source_contents, processing_contexts)

    emit:
    rcv_xml = clinvar_rcv
    vcv_xml = clinvar_vcv
    records = PARSE_CLINVAR.out.records
    summary = PARSE_CLINVAR.out.summary
    warnings = PARSE_CLINVAR.out.warnings
    rejections = PARSE_CLINVAR.out.rejections
}
