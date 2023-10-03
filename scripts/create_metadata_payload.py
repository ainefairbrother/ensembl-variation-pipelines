#!/usr/bin/env python3

import sys
import configparser
import argparse
import os
import json
import subprocess
import requests

ENDPOINT = "https://services.test.ensembl-production.ebi.ac.uk/api/genome_metadata/datasets/"

def parse_args(args = None):
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--api_outdir", dest="api_outdir", type=str, help="path to a vcf prepper api output directory")
    parser.add_argument("--endpoint", dest="endpoint", type=str, help="metadata api url")
    parser.add_argument("--debug", dest="debug", action="store_true")
    
    return parser.parse_args(args)
    
def get_variant_count(file: str) -> str:
    process = subprocess.run(["bcftools", "index", "--nrecords", file],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )

    try:
        return int(process.stdout.decode().strip())
    except Exception as e:
        print(f"""Could not get count from {file}
        {e}""")
        return None

def submit_payload(endpoint: str, payload: str) -> str:
    requests.put(endpoint, payload)
    
def main(args = None):
    args = parse_args(args)
    
    api_outdir = args.api_outdir or os.getcwd()
    endpoint = args.endpoint or ENDPOINT
    debug = args.debug

    for genome_uuid in os.listdir(api_outdir):
        api_vcf = os.path.join(api_outdir, genome_uuid, "variation.vcf.gz")
        if not os.path.isfile(api_vcf):
            continue

        variant_count = get_variant_count(api_vcf)
        
        if variant_count is not None:
            payload = {}
            payload["user"] = "nakib"
            payload["name"] = "variation"
            payload["description"] = "A test variation dataset",
            payload["label"] = "Manual Add",
            payload["dataset_type"] = "variation", 
            
            dataset_source = {}
            dataset_source["name"] = api_vcf
            dataset_source["type"] = "vcf"
            payload["dataset_source"] = dataset_source
            
            payload["genome_uuid"] = genome_uuid

            dataset_attribute = {}
            dataset_attribute["name"] = "short_variants"
            dataset_attribute["value"] = variant_count
            payload["dataset_attribute"] = dataset_attribute
            
            if not debug:
                submit_payload(endpoint, payload)
            else:
                print(json.dumps(payload, indent = 4))
    
if __name__ == "__main__":
    sys.exit(main())