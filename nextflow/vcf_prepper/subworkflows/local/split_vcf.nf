//
// Split large VCF
//

include { READ_VCF_CHR } from "../../modules/local/read_vcf_chr.nf"
include { SPLIT_VCF_CHR } from "../../modules/local/split_vcf_chr.nf"

workflow SPLIT_VCF {
  take:
    input
  
  main:
    READ_VCF_CHR( input )
    SPLIT_VCF_CHR(READ_VCF_CHR.out.transpose())

  emit:
    SPLIT_VCF_CHR.out
}