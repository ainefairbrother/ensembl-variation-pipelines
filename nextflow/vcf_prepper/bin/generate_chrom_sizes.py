#!/usr/bin/env python3

import sys
import configparser
import argparse
import subprocess
import os

from helper import parse_ini, get_db_name

def parse_args(args = None):
    parser = argparse.ArgumentParser()
    
    parser.add_argument(dest="genome", type=str, help="genome string with format <species>_<assembly>")
    parser.add_argument(dest="version", type=int, help="Ensembl release version")
    parser.add_argument('-I', '--ini_file', dest="ini_file", type=str, required = False, help="full path database configuration file, default - DEFAULT.ini in the same directory.")
    parser.add_argument('--chrom_sizes', dest="chrom_sizes", type=str, required = False, help="file with chromomsome sizes, default - <genome>.chrom.sizes in the same directory.")
    parser.add_argument('--force', dest="force", action="store_true", help="forcefully create config even if already exists")
    
    return parser.parse_args(args)

def generate_chrom_sizes(server: dict, core_db: str, chrom_sizes: str, assembly: str = "grch38", force: bool = False) -> None:
    if os.path.exists(chrom_sizes) and not force:
        print(f"[INFO] {chrom_sizes} file already exists, skipping ...")
        return
        
    query = f"SELECT coord_system_id FROM coord_system WHERE version = '{assembly}';"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "--database", core_db,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    coord_ids = tuple([int(id) for id in process.stdout.decode().strip().split("\n")])
    
    query = f"SELECT name, length FROM seq_region WHERE coord_system_id IN {coord_ids};"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "--database", core_db,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    with open(chrom_sizes, "w") as file: 
        file.write(process.stdout.decode())
        
    query = f"SELECT ss.synonym, s.length FROM seq_region AS s, seq_region_synonym AS ss WHERE s.seq_region_id = ss.seq_region_id AND s.coord_system_id IN {coord_ids};"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "--database", core_db,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    with open(chrom_sizes, "a") as file: 
        file.write(process.stdout.decode().strip())
    
    # remove duplicates    
    with open(chrom_sizes, "r") as file: 
        lines = file.readlines()
        
    lengths = {}
    for line in lines:
        name, length = [col.strip() for col in line.split("\t")]  
        if name not in lengths or int(lengths[name]) < int(length): 
            lengths[name] = length
            
    with open(chrom_sizes, "w") as file:
        for name in lengths:
            # we will keep length + 1 because bedToBigBed fails if it finds variant at boundary
            length = int(lengths[name]) + 1
            file.write(f"{name}\t{str(length)}\n")

def main(args = None):
    args = parse_args(args)
    
    chrom_sizes = args.chrom_sizes or f"{args.genome}.chrom.sizes"
    assembly = args.genome.split("_")[-1]
    species = args.genome.replace(f"_{assembly}", "")
    db_server = parse_ini(args.ini_file, assembly)
    core_db = get_db_name(db_server, args.version, species, type = "core")
    
    generate_chrom_sizes(db_server, core_db, chrom_sizes, assembly, args.force)
    
if __name__ == "__main__":
    sys.exit(main())