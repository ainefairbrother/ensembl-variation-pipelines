#!/usr/bin/env python3

import pytest
import os
import sys
import json
from uuid import UUID
import argparse
import configparser
import subprocess
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def parse_args(args = None):
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--dir", dest="dir", type=str, default = os.getcwd())
    parser.add_argument("--input_config", dest="input_config", type=str)
    parser.add_argument("-O", "--output_dir", dest="output_dir", type=str, default = os.getcwd())
    
    return parser.parse_args(args)

def get_species_metadata(input_config: str = None) -> dict:
    if input_config is None or not os.path.isfile(input_config):
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
            species_metadata[genome_uuid]["file_location"] = genome["file_location"]

    return species_metadata

def is_valid_uuid(uuid: str):
    try:
        uuid_obj = UUID(uuid)
    except ValueError:
        return False
    return str(uuid_obj) == uuid

def main(args = None):
    args = parse_args(args)

    input_config = args.input_config or None
    dir = args.dir
    output_dir = args.output_dir

    # create output directory if not already exist
    os.makedirs(output_dir, exist_ok=True)

    species_metadata = {}
    if input_config is not None:
        species_metadata = get_species_metadata(input_config)

    vcf_files = []
    api_outdir = os.path.join(dir, "api")
    track_outdir = os.path.join(dir, "tracks")
    for genome_uuid in os.listdir(api_outdir):
        if species_metadata and genome_uuid not in species_metadata:
            continue

        if not is_valid_uuid(genome_uuid):
            logger.warning(f"{genome_uuid} is not a valid uuid")
            continue

        species = species_metadata[genome_uuid]["species"]
        vcf = os.path.join(api_outdir, genome_uuid, "variation.vcf.gz")
        source_vcf = species_metadata[genome_uuid]["file_location"]
        bigbed = os.path.join(track_outdir, genome_uuid, "variant-details.bb")
        bigwig = os.path.join(track_outdir, genome_uuid, "variant-summary.bw")

        subprocess.run([
                "bsub",
                "-J", f"dc_{species}",
                "-oo", f"{output_dir}/dc_{species}.out",
                "-eo", f"{output_dir}/dc_{species}.err",
                f"pytest " + \
                f"--source_vcf={source_vcf} " + \
                f"--bigbed={bigbed} " + \
                f"--bigwig={bigwig} " + \
                f"--vcf={vcf} " + \
                f"--species={species} " + \
                "./"
            ]
        )

if __name__ == "__main__":
    sys.exit(main())