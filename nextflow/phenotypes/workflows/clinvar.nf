include { DOWNLOAD_CLINVAR } from '../modules/local/download_clinvar'
include { PARSE_CLINVAR } from '../modules/local/parse_clinvar'

workflow CLINVAR {
    main:
    releases_to_download = []
    if (!params.rcv_input) {
        releases_to_download.add('RCV')
    }
    if (!params.vcv_input) {
        releases_to_download.add('VCV')
    }
    if (releases_to_download) {
        DOWNLOAD_CLINVAR(channel.fromList(releases_to_download))
    }

    clinvar_rcv = params.rcv_input ? channel.fromPath(params.rcv_input, checkIfExists: true) :
        DOWNLOAD_CLINVAR.out.xml.filter { release_type, xml -> release_type == 'RCV' }.map { release_type, xml -> xml }
    clinvar_vcv = params.vcv_input ? channel.fromPath(params.vcv_input, checkIfExists: true) :
        DOWNLOAD_CLINVAR.out.xml.filter { release_type, xml -> release_type == 'VCV' }.map { release_type, xml -> xml }

    // Stage the package so edits invalidate parsing, not completed downloads.
    python_source = file("${projectDir}/src", checkIfExists: true)
    PARSE_CLINVAR(clinvar_rcv, clinvar_vcv, python_source, params.assembly)

    emit:
    rcv_xml = clinvar_rcv
    vcv_xml = clinvar_vcv
    records = PARSE_CLINVAR.out.records
    summary = PARSE_CLINVAR.out.summary
    warnings = PARSE_CLINVAR.out.warnings
    rejections = PARSE_CLINVAR.out.rejections
}
