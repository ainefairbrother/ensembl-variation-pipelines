#!/usr/bin/env python3

import sys
from cyvcf2 import VCF, Writer
import argparse

META = """##fileformat=VCFv4.2
##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of the variation data">
"""
HEADER="""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def parse_args(args = None, description: bool = None):
    parser = argparse.ArgumentParser(description = description)
    
    parser.add_argument(dest="input_file", type=str, help="input VCF file")
    parser.add_argument(dest="source", type=str, help="input VCF file source")
    parser.add_argument(dest="synonym_file", type=str, help="text file with chrmosome synonyms")
    parser.add_argument('--rename_clinvar_ids', dest="rename_clinvar_ids", action="store_true")
    parser.add_argument('--chromosomes', dest="chromosomes", type=str, help="comma separated list of chromosomes to put in header")
    parser.add_argument('-O', '--output_file', dest="output_file", type=str)
    
    return parser.parse_args(args)

def format_clinvar_id(id: str) -> str:
    if not id.startswith("VCV"):
        leading_zero = 9 - len(id)
        return "VCV" + ("0" * leading_zero) + id

    return id

def format_meta(meta: str, chromosomes: str = None, synonyms: list = None) -> str:
    if chromosomes is None:
        return meta

    for chromosome in chromosomes.split(","):
        chr_syn = synonyms[chromosome] if chromosome in synonyms else chromosome
        meta += f"##contig=<ID={chr_syn}>\n"
    return meta

def main(args = None):
    args = parse_args(args)

    input_file = args.input_file
    source = args.source
    synonym_file = args.synonym_file
    chromosomes = args.chromosomes or None
    output_file = args.output_file or input_file.replace(".vcf.gz", "_renamed.vcf.gz")

    synonyms = {}
    with open(synonym_file) as file:
        for line in file:
            chr = line.split("\t")[0].strip() 
            synonym = line.split("\t")[1].strip()

            synonyms[chr] = synonym

    if args.rename_clinvar_ids and source == "ClinVar":
        format_id = format_clinvar_id
    else:
        format_id = lambda x : x
    
    meta = format_meta(META, chromosomes, synonyms)

    with open(output_file, "w") as o_file:
        o_file.write(meta)
        o_file.write(HEADER)

        input_vcf = VCF(input_file)
        for variant in input_vcf:
            o_file.write("\t".join([
                    synonyms[variant.CHROM] if variant.CHROM in synonyms else variant.CHROM,
                    str(variant.POS),
                    format_id(variant.ID),
                    variant.REF,
                    ",".join(variant.ALT),
                    ".",
                    ".",
                    f"SOURCE={source}"
                ]) + "\n"
            )
        input_vcf.close()
        
    o_file.close()
    
if __name__ == "__main__":
    sys.exit(main())