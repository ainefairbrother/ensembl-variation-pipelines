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
    parser.add_argument("--input_config", dest="input_config", type=str, help="input_config json file used in vcf_prepper")
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

def parse_input_config(input_config: str) -> dict:
    if not os.path.isfile(input_config):
        return []

    with open(input_config, "r") as file:
        input_config_json = json.load(file)

    species_metadata = {}
    for species in input_config_json:
        for genome in input_config_json[species]:
            genome_uuid = genome["genome_uuid"]
            if genome_uuid not in species_metadata:
                species_metadata[genome_uuid] = {}

            species_metadata[genome_uuid]["species"] = genome["species"]
            species_metadata[genome_uuid]["assembly"] = genome["assembly"]

    return species_metadata

def submit_payload(endpoint: str, payload: str) -> str:
    requests.put(endpoint, payload)
    
def main(args = None):
    args = parse_args(args)
    
    api_outdir = args.api_outdir or os.getcwd()
    input_config = args.input_config or None
    endpoint = args.endpoint or None
    debug = args.debug

    if not debug and endpoint is None:
        print("[ERROR] please provide an endpoint using --endpoint if not using debug mode")
        exit(1)

    species_metadata = {}
    if input_config is not None:
        species_metadata = parse_input_config(input_config)

    print(f"[INFO] checking directory - {api_outdir} for genome uuids")
    for genome_uuid in os.listdir(api_outdir):
        if species_metadata and genome_uuid not in species_metadata:
            continue

        if not is_valid_uuid(genome_uuid):
            print(f"[WARN] {genome_uuid} is not a valid uuid")
            continue

        api_vcf = os.path.join(api_outdir, genome_uuid, "variation.vcf.gz")
        if not os.path.isfile(api_vcf):
            print(f"[WARN] file not found - {api_vcf}")
            continue

        # TBD: get this data from thoas if input_config not given
        species = species_metadata[genome_uuid]["species"]
        assembly = species_metadata[genome_uuid]["assembly"]
        variant_count = get_variant_count(api_vcf)
        
        if variant_count is not None:
            payload = {}
            payload["user"] = "nakib"
            payload["name"] = "variation"
            payload["description"] = f"Short variant data for {species}"
            payload["label"] = assembly
            payload["dataset_type"] = "variation"
            
            dataset_source = {}
            dataset_source["name"] = api_vcf
            dataset_source["type"] = "vcf"
            payload["dataset_source"] = dataset_source
            
            payload["genome_uuid"] = genome_uuid

            dataset_attribute = {}
            dataset_attribute["name"] = "variation.short_variants"
            dataset_attribute["value"] = variant_count
            payload["dataset_attribute"] = dataset_attribute
            
            if not debug:
                submit_payload(endpoint, payload)
            else:
                print(json.dumps(payload, indent = 4))
    
if __name__ == "__main__":
    sys.exit(main())