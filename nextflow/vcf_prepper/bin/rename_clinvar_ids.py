#!/usr/bin/env python3

import sys
from cyvcf2 import VCF, Writer
import argparse

def parse_args(args = None, description: bool = None):
    parser = argparse.ArgumentParser(description = description)
    
    parser.add_argument(dest="input_file", type=str, help="input VCF file")
    parser.add_argument('-O', '--output_file', dest="output_file", type=str)
    
    return parser.parse_args(args)
    
def get_accession_id(id: str) -> str:
    leading_zero = 9 - len(id)
    return "VCV" + ("0" * leading_zero) + id

def main(args = None):
    args = parse_args(args)

    input_file = args.input_file
    output_file = args.output_file or input_file.replace(".vcf.gz", "_renamed.vcf.gz")
    
    output_vcf_writer = Writer(ouput_file, input_vcf)
    for variant in input_vcf:
        id = variant.ID
        
        if not id.startswith("VCV"):
            variant.ID = get_accession_id(id)       
        
        output_vcf_writer.write_record(variant)
    output_vcf_writer.close()
    input_vcf.close()
    
if __name__ == "__main__":
    sys.exit(main())