#!/usr/bin/env python3

import sys
import configparser
import argparse
import os
import json
import subprocess
import requests
from uuid import UUID
from cyvcf2 import VCF

def parse_args(args = None):
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--input_config", dest="input_config", type=str, required = True, help="input_config json file used in vcf_prepper")
    
    return parser.parse_args(args)

def is_valid_uuid(uuid: str):
    try:
        uuid_obj = UUID(uuid)
    except ValueError:
        return False
    return str(uuid_obj) == uuid

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

            species_metadata[genome_uuid]["source_name"] = genome["source_name"]
            species_metadata[genome_uuid]["species"] = genome["species"]

    return species_metadata

def get_source_info(source: str) -> str:
    if source == "dbSNP":
        return "dbSNP - build 156"
    elif source == "EVA":
        return "European Variation Archive (EVA) - release 5"
    elif source == "Ensembl":
        return "Ensembl - e110"

def get_source_url(source: str) -> str:
    if source == "dbSNP":
        return "https://www.ncbi.nlm.nih.gov/snp"
    elif source == "EVA":
        return "https://www.ebi.ac.uk/eva"
    elif source == "Ensembl":
        return "https://www.ensembl.org/index.html"
    
def main(args = None):
    args = parse_args(args)
    
    input_config = args.input_config

    species_metadata = {}
    if input_config is not None:
        species_metadata = parse_input_config(input_config)

    metadata = {}
    for genome_uuid in species_metadata:
        if not is_valid_uuid(genome_uuid):
            print(f"[WARN] {genome_uuid} is not a valid uuid")
            continue

        source = species_metadata[genome_uuid]["source_name"]
        species = species_metadata[genome_uuid]["species"]

        source_info = get_source_info(source)
        source_url = get_source_url(source)

        metadata[genome_uuid] = {}
        metadata[genome_uuid]["label"] = f"{source} short variants"

        metadata[genome_uuid]["datafiles"] = {}
        metadata[genome_uuid]["datafiles"]["details"] = f"variant-{source.lower()}-details.bb"
        metadata[genome_uuid]["datafiles"]["summary"] = f"variant-{source.lower()}-summary.bw"

        metadata[genome_uuid]["description"] = f"All short variants (SNPs and indel) data from {source_info}"
        
        metadata[genome_uuid]["source"] = {}
        metadata[genome_uuid]["source"]["name"] = source
        metadata[genome_uuid]["source"]["url"] = source_url 
            
    print(json.dumps(metadata, indent = 4))
    
if __name__ == "__main__":
    sys.exit(main())