include { PHENOTYPES } from './workflows/phenotypes'

workflow {
    if (params.containsKey('input') || params.containsKey('rcv_input') || params.containsKey('vcv_input')) {
        error('Input files belong in sources.<source>.inputs in --source_config, not pipeline CLI options. Omitted inputs are downloaded.')
    }
    if (params.containsKey('assembly') || params.containsKey('env_file')) {
        error('Use --source_config for species and assembly settings. Parsing no longer queries metadata or reads .env.')
    }
    PHENOTYPES()
}
