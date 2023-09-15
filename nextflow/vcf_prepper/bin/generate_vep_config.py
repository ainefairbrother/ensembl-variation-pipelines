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
    parser.add_argument(dest="version", type=str, help="Ensembl release version")
    parser.add_argument('-I', '--ini_file', dest="ini_file", type=str, required = False, help="full path database configuration file, default - DEFAULT.ini in the same directory.")
    parser.add_argument('--vep_config', dest="vep_config", type=str, required = False, help="vep configuration file, default - <genome>.ini in the same directory.")
    parser.add_argument('--force', dest="force", action="store_true", help="forcefully create config even if already exists")
    
    return parser.parse_args(args)

def generate_vep_config(server: dict, db_name: str, vep_config: str, force: bool = False) -> None:
    if os.path.exists(synonym_file) and not force:
        print(f"[INFO] {synonym_file} file already exists, skipping ...")
        return
        
    query = f"SELECT ss.synonym, sr.name FROM seq_region AS sr, seq_region_synonym AS ss WHERE sr.seq_region_id = ss.seq_region_id;"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "--database", db_name,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    with open(synonym_file, "w") as file: 
        file.write(process.stdout.decode().strip())
    
    # remove duplicates and change seq region name that are longer than 31 character    
    with open(synonym_file, "r") as file: 
        lines = file.readlines()
        
    names = {}
    for line in lines:
        synonym, name = [col.strip() for col in line.split("\t")]  
        if synonym not in names or len(names[synonym]) > len(name): 
            names[synonym] = name
    
    # if name is longer than 31 character we try to take a synonym
    new_names = {}
    for synonym in names:
        name = names[synonym]
        if len(name) > 31:
            # if the current synonym is less than 31 character we do not need to have it in the file
            if len(synonym) <= 31:
                pass
            # if the current synonym is longer than 31 character we look for other synonym of the name
            else:
                change_name = synonym
                for alt_synonym in names:
                    if names[alt_synonym] == synonym and len(alt_synonym) < 31:
                        change_name = alt_synonym
                        
                if len(change_name) > 31:
                    print(f"[WARNING] cannot resolve {name} to a synonym which is under 31 character")
                            
                new_names[synonym] = change_name        
        else:
            new_names[synonym] = name
    
    
    with open(synonym_file, "w") as file:
        for synonym in new_names:
            file.write(f"{synonym}\t{new_names[synonym]}\n")
    
def main(args = None):
    args = parse_args(args)
    
    vep_config = args.vep_config or f"{args.genome}.txt"
    assembly = args.genome.split("_")[-1]
    species = args.genome.replace(f"_{assembly}", "")

    db_server = parse_ini(args.ini_file, assembly)
    db_name = get_db_name(db_server, args.version, species)
    
    generate_vep_config(db_server, db_name, vep_config, args.force)
    
if __name__ == "__main__":
    sys.exit(main())