#!/usr/bin/env python3

import sys
import configparser
import argparse
import subprocess
import os
import requests
import glob

from helper import *

FASTA_DIR = "/nfs/production/flicek/ensembl/variation/data/VEP/fasta"

def parse_args(args = None):
    parser = argparse.ArgumentParser()
    
    parser.add_argument(dest="species", type=str, help="species production name")
    parser.add_argument(dest="assembly", type=str, help="assembly default")
    parser.add_argument(dest="version", type=int, help="Ensembl release version")
    parser.add_argument('--division', dest="division", type=str, required = False, help="Ensembl division the species belongs to")
    parser.add_argument('-I', '--ini_file', dest="ini_file", type=str, required = False, help="full path database configuration file, default - DEFAULT.ini in the same directory.")
    parser.add_argument('--fasta_dir', dest="fasta_dir", type=str, required = False, help="FASTA directory")
    
    return parser.parse_args(args)
    
def ungzip_fasta(fasta_dir: str, compressed_fasta: str) -> str:
    if os.path.dirname(compressed_fasta) != fasta_dir:
        print(f"[ERROR] Fasta file {fasta_dir} in wrong directory; should be in - {fasta_dir}")
        exit(1)
        
    process = subprocess.run(["gzip", "-d", compressed_fasta],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    if process.returncode != 0:
        print(f"[ERROR] Could not uncompress fasta file - {compressed_fasta}")
        exit(1)
        
    return compressed_fasta[:-3]
    
def bgzip_fasta(fasta_dir: str, unzipped_fasta: str) -> None:
    if os.path.dirname(unzipped_fasta) != fasta_dir:
        print(f"[ERROR] Fasta file {fasta_dir} in wrong directory; should be in - {fasta_dir}")
        exit(1)
        
    process = subprocess.run(["bgzip", unzipped_fasta],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    if process.returncode != 0:
        print(f"[ERROR] Could not bgzip fasta file - {unzipped_fasta}")
        exit(1)
    
def main(args = None):
    args = parse_args(args)
    
    species = args.species
    assembly = args.assembly
    version = args.version
    ini_file = args.ini_file or "DEFAULT.ini"
    core_server = parse_ini(ini_file, "core")
    core_db = get_db_name(core_server, args.version, species, type = "core")
    division = args.division or get_division(core_db, core_db)
    species_url_name = get_species_url_name(core_db, core_db)
    
    fasta_dir = args.fasta_dir or FASTA_DIR
    fasta_glob = os.path.join(fasta_dir, f"{species_url_name}.{assembly}.dna.*.fa.gz")
    if not glob.glob(fasta_glob):
        print(f"[INFO] {fasta_glob} does not exists. Creating ...")
        
        rl_version = get_relative_version(version, division)
        src_compressed_fasta = get_ftp_path(species, assembly, division, rl_version, "fasta", "local", species_url_name)
        
        if src_compressed_fasta is not None:
            compressed_fasta = os.path.join(fasta_dir, os.path.basename(src_compressed_fasta))
            returncode = copyto(src_compressed_fasta, compressed_fasta)
        
        if src_compressed_fasta is None or returncode != 0:
            print(f"[INFO] Failed to copy fasta file - {src_compressed_fasta}, will retry using remote FTP")
            
            compressed_fasta_url = get_ftp_path(species, assembly, division, rl_version, "fasta", "remote", species_url_name)
            
            compressed_fasta = os.path.join(fasta_dir, compressed_fasta_url.split('/')[-1])
            returncode = download_file(compressed_fasta, compressed_fasta_url)
            if returncode != 0:
                print(f"[ERROR] Could not download fasta file - {compressed_fasta_url}")
                exit(1)
        
        unzipped_fasta = ungzip_fasta(fasta_dir, compressed_fasta)
        bgzip_fasta(fasta_dir, unzipped_fasta)
        
if __name__ == "__main__":
    sys.exit(main())