#!/usr/bin/env python3

import subprocess
import configparser
import os
import pwd

def parse_ini(ini_file: str, species: str = "homo_sapiens", assembly: str = "grch38") -> dict:
    config = configparser.ConfigParser()
    config.read(ini_file)
    
    if species == "homo_sapiens" and assembly.lower() == "grch37":
        if not "grch37_database" in config:
            print(f"[ERROR] Could not find 'database' config in ini file - {ini_file}")
            exit(1)
        else:
            host = config["grch37_database"]["host"]
            port = config["grch37_database"]["port"]
            user = config["grch37_database"]["user"]
    else:
        if not "database" in config:
            print(f"[ERROR] Could not find 'database' config in ini file - {ini_file}")
            exit(1)
        else:
            host = config["database"]["host"]
            port = config["database"]["port"]
            user = config["database"]["user"]
    
    return {
        "host": host, 
        "port": port, 
        "user": user
    }

def get_db_name(server: dict, version: str, species: str = "homo_sapiens", type: str = "core") -> str:
    query = f"SHOW DATABASES LIKE '{species}_{type}%{version}%';"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    return process.stdout.decode().strip()

def get_division(server: dict, core_db: str) -> str:
    query = "SELECT meta_value FROM meta WHERE meta_key = 'species.division';"
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
    return process.stdout.decode().strip()

def get_species_url_name(server: dict, core_db: str) -> str:
    query = "SELECT meta_value FROM meta WHERE meta_key = 'species.url';"
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
    return process.stdout.decode().strip()
    
def get_relative_version(version: int, division: str = "EnsemblVertebrates") -> int:
    return (version - 53) if division != "EnsemblVertebrates" else version
    
def download_file(local_filename: str, url: str) -> int:
    process = subprocess.run(["wget", url, "-O", local_filename],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
        
    return process.returncode 
    
def get_ftp_path(
        species: str, 
        assembly: str, 
        division: str, 
        version: int, 
        type: str = "cache", 
        mode: str = "local",
        species_url_name: str = None
    ) -> str:
    version = str(version)
    
    base = "/nfs/production/flicek/ensembl/production/ensemblftp" if mode == "local" else "ftp.ensembl.org/pub"
    release_segment = f"release-{version}"
    
    division_segment = ""
    if division != "EnsemblVertebrates" and type != "conservation":
        division_segment = f"{division[7:].lower()}"
        
    if type == "cache":
        prefix = "variation/indexed_vep_cache"
    elif type == "fasta":
        prefix = f"fasta/{species}/dna"
    elif type == "conservation":
        prefix = "compara/conservation_scores/91_mammals.gerp_conservation_score"
    
    if type == "cache":
        file_name = f"{species}_vep_{version}_{assembly}.tar.gz"
    elif type == "fasta":
        file_name = f"{species_url_name}.{assembly}.dna.toplevel.fa.gz"
    elif type == "conservation":
        file_name = f"gerp_conservation_scores.{species}.{assembly}.bw"

    full_path = os.path.join(base, release_segment, division_segment, prefix, file_name)
    
    if mode == "local" and os.path.isfile(full_path):
        return full_path
    elif mode == "remote":
        return f"https://{full_path}"
    
    print(full_path)
    return None
    
def copyto(src_file: str, dest_file: str) -> int:
    process = subprocess.run(["rsync", src_file, dest_file], 
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
        
    return process.returncode