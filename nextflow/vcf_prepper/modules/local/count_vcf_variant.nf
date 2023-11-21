#!/usr/bin/env nextflow

process COUNT_VCF_VARIANT {
    label 'bcftools'

    input:
    tuple val(meta), val(vcf), val(vcf_index)

    output:
    tuple val(meta), val(vcf), val(vcf_index), env(count)

    shell:
    '''
    count=$(bcftools index --nrecords !{vcf})
    '''
}