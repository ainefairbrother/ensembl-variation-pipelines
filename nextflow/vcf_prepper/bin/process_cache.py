#!/usr/bin/env python3

import sys
import configparser
import argparse
import subprocess
import os
import requests

from helper import *

CACHE_DIR = "/nfs/production/flicek/ensembl/variation/data/VEP/tabixconverted"

def parse_args(args = None):
    parser = argparse.ArgumentParser()
    
    parser.add_argument(dest="genome", type=str, help="genome string with format <species>_<assembly>")
    parser.add_argument(dest="version", type=int, help="Ensembl release version")
    parser.add_argument('--division', dest="division", type=str, required = False, help="Ensembl division the species belongs to")
    parser.add_argument('-I', '--ini_file', dest="ini_file", type=str, required = False, help="full path database configuration file, default - DEFAULT.ini in the same directory.")
    parser.add_argument('--cache_dir', dest="cache_dir", type=str, required = False, help="VEP cache directory")
    
    return parser.parse_args(args)
    
def uncompress_cache(cache_dir: str, compressed_cache: str) -> None:
    process = subprocess.run(["tar", "-xvzf", compressed_cache, "-C", cache_dir],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    if process.returncode != 0:
        print(f"[ERROR] Could not uncompress cache file - {compressed_cache}")
        exit(1)
    
def main(args = None):
    args = parse_args(args)
    
    version = args.version
    assembly = args.genome.split("_")[-1]
    species = args.genome.replace(f"_{assembly}", "")
    db_server = parse_ini(args.ini_file, assembly)
    core_db = get_db_name(db_server, args.version, species, type = "core")
    division = args.division or get_division(db_server, core_db)
    
    cache_dir = args.cache_dir or CACHE_DIR
    rl_version = get_relative_version(version, division)
    genome_cache_dir = os.path.join(cache_dir, species, f"{rl_version}_{assembly}")         
    if not os.path.exists(genome_cache_dir):
        print(f"[INFO] {genome_cache_dir} directory does not exists. Creating ...")
        
        compressed_cache = get_ftp_path(species, assembly, division, rl_version, "cache")
        
        if compressed_cache is None:
            print(f"[INFO] Could not find cache in local ftp directory, will retry using remote FTP")
            
            compressed_cache_url = get_ftp_path(species, assembly, division, rl_version, "cache", "remote")
            
            compressed_cache = os.path.join(cache_dir, compressed_cache_url.split('/')[-1])
            returncode = download_file(compressed_cache, compressed_cache_url)
            if returncode != 0:
                print(f"[ERROR] Could not download cache file - {compressed_cache_url}")
                exit(1)
        
        uncompress_cache(cache_dir, compressed_cache)   
    
if __name__ == "__main__":
    sys.exit(main())