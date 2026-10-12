include { CLINVAR } from './clinvar'

def read_source_config(config_path) {
    def config_file = file(config_path, checkIfExists: true).toAbsolutePath()
    def config = new groovy.json.JsonSlurper().parseText(config_file.text)
    if (!(config instanceof Map) || !(config.sources instanceof Map) || !config.sources) {
        error('Source configuration must contain a non-empty sources object.')
    }
    config.sources.each { source, settings ->
        // Add a workflow here when another source is implemented; never silently skip it.
        if (source != 'clinvar') {
            error("No workflow is implemented for source '${source}'.")
        }
        if (!(settings instanceof Map) || !(settings.contexts instanceof List) || !settings.contexts) {
            error("Source '${source}' must contain a non-empty contexts list.")
        }
        def inputs = settings.containsKey('inputs') ? settings.inputs : [:]
        if (!(inputs instanceof Map)) {
            error("Source '${source}' inputs must be an object.")
        }
        settings.inputs = inputs.collectEntries { name, input_path ->
            if (!(name in ['rcv', 'vcv'])) {
                error("Unknown ClinVar input '${name}'; use 'rcv' or 'vcv'.")
            }
            if (input_path != null && (!(input_path instanceof String) || !input_path.trim())) {
                error("sources.${source}.inputs.${name} must be a file path or null.")
            }
            // Absolute paths stay absolute; relative paths belong to the
            // configuration directory, not the directory Nextflow launched in.
            return [name, input_path == null ? null : config_file.parent.resolve(input_path).normalize().toString()]
        }
        def seen = [] as Set
        settings.contexts.each { context ->
            if (!(context instanceof Map)
                || !(context.species_taxonomy_id instanceof Integer)
                || context.species_taxonomy_id != 9606
                || context.species_production_name != 'homo_sapiens'
                || !(context.assembly instanceof String)
                || !(context.assembly ==~ /[A-Za-z0-9_.-]+/)
                || !(context.assembly_uuid instanceof String)
                || !(context.assembly_uuid ==~ /[0-9a-fA-F]{8}(-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}/)) {
                error('Each ClinVar context requires species_taxonomy_id 9606, species_production_name homo_sapiens, an assembly label and an assembly UUID.')
            }
            def key = "${context.species_production_name}/${context.assembly}".toString()
            if (!seen.add(key)) {
                error("Duplicate context '${source}/${key}' in source configuration.")
            }
        }
    }
    return config
}

workflow PHENOTYPES {
    main:
    source_config = read_source_config(params.source_config)
    CLINVAR(source_config.sources.clinvar)

    emit:
    records = CLINVAR.out.records
    summary = CLINVAR.out.summary
    warnings = CLINVAR.out.warnings
    rejections = CLINVAR.out.rejections
}
