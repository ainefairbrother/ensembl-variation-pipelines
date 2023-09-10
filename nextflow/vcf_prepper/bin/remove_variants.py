#!/usr/bin/env python3

import sys
from cyvcf2 import VCF, Writer
from cyvcf2.cyvcf2 import Variant
import argparse
from argparse import RawTextHelpFormatter
from typing import Callable

def parse_args(args = None, description: bool = None):
    parser = argparse.ArgumentParser(description = description, formatter_class=RawTextHelpFormatter)
    
    parser.add_argument(dest="input_file", type=str, help="input VCF file")
    parser.add_argument('--remove_nonunique_ids', dest="remove_nonunique_ids", action="store_true", help="remove variants with same ids")
    parser.add_argument('--remove_patch_regions', dest="remove_patch_regions", action="store_true", help="remove variant in patch region")
    parser.add_argument('-O', '--output_file', dest="output_file", type=str)
    
    return parser.parse_args(args)
 
def get_id(variant: Variant) -> str:
    'Get variant id'
    
    return variant.ID

def get_positioned_id(variant: Variant) -> str:
    'Get variant positioned id'
        
    id = variant.ID or "unknown"
    return variant.CHROM + ":" + str(variant.POS) + ":" + id

def generate_removal_status(vcf_file: str, get_identifier: Callable, remove_patch_regions: bool = True) -> dict:
    'Generate hash against variant about its removal status'
    
    removal_status = {}
    input_vcf = VCF(vcf_file)
    for variant in input_vcf:
        variant_identifier = get_identifier(variant)
        # Order is important here. Check for uniqueness is based on existance - we should check it first
        removal_status[variant_identifier] = variant_identifier in removal_status
        if remove_patch_regions:
            chr = variant.CHROM
            removal_status[variant_identifier] = removal_status[variant_identifier] or ("CTG" in chr) or ("PATCH" in chr) or ("TEST" in chr)
    input_vcf.close()
    
    return removal_status

def main(args = None):
    description = '''
    Removes variant based on uniqueness and sequence region. 
        1) By default, variant is discarded if the positioned identifier (chrom:position:id) is same for multiple variant record. The assumption is that the variants will be multi-allelic if needed be instead of bi-allelic in the source VCF file.
        2) Optionally, we can ask to remove variant with same ids even if they are in different location (using the remove_nonunique_ids argument).
        3) When removed, all the variant record is removed. For example, if there is two variant record with same positioned id then both of them will be removed.
    '''
    args = parse_args(args, description)
    
    input_file = args.input_file
    remove_nonunique_ids = args.remove_nonunique_ids or False
    remove_patch_regions = args.remove_patch_regions or True
    output_file = args.output_file or input_file.replace("renamed", "processed")
    
    if remove_nonunique_ids:
        get_identifier = get_id
    else:
        get_identifier = get_positioned_id
        
    removal_status = generate_removal_status(input_file, get_identifier, remove_patch_regions)
    
    # Remove variant based on removal status
    input_vcf = VCF(input_file)
    output_vcf_writer = Writer(output_file, input_vcf)
    for variant in input_vcf:
        
        variant_identifier = get_identifier(variant)
        if removal_status[variant_identifier]:
                continue
                
        output_vcf_writer.write_record(variant)
    output_vcf_writer.close()
    input_vcf.close()
    
if __name__ == "__main__":
    sys.exit(main())