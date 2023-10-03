#!/usr/bin/env python3

import sys
import configparser
import argparse
import os
import json
import subprocess
import requests
from uuid import UUID

def parse_args(args = None):
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--api_outdir", dest="api_outdir", type=str, help="path to a vcf prepper api output directory")
    parser.add_argument("--endpoint", dest="endpoint", type=str, help="metadata api url")
    parser.add_argument("--debug", dest="debug", action="store_true")
    
    return parser.parse_args(args)

def is_valid_uuid(uuid: str):
    try:
        uuid_obj = UUID(uuid)
    except ValueError:
        return False
    return str(uuid_obj) == uuid

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
    endpoint = args.endpoint or None
    debug = args.debug

    if not debug and endpoint is None:
        print("[ERROR] please provide an endpoint using --endpoint if not using debug mode")
        exit(1)

    print(f"[INFO] checking directory - {api_outdir} for genome uuids")
    for genome_uuid in os.listdir(api_outdir):
        if not is_valid_uuid(genome_uuid):
            print(f"[WARN] {genome_uuid} is not a valid uuid")
            continue

        api_vcf = os.path.join(api_outdir, genome_uuid, "variation.vcf.gz")
        if not os.path.isfile(api_vcf):
            print(f"[WARN] file not found - {api_vcf}")
            continue

        variant_count = get_variant_count(api_vcf)
        
        if variant_count is not None:
            payload = {}
            payload["user"] = "nakib"
            payload["name"] = "variation"
            payload["description"] = "Dataset with short variant data for TBD" # add species name when thoas is up 
            payload["label"] = "TBD"        # add assembly accession when thoas is up
            payload["dataset_type"] = "variation"
            
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