#!/usr/bin/env nextflow

process COUNT_VCF_VARIANT {
    label 'bcftools'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path(vcf), path(vcf_index), env(count)

    shell:
    '''
    count=$(bcftools index --nrecords !{vcf})
    '''
}