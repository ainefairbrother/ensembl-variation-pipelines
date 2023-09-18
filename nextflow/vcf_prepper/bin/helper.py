#!/usr/bin/env python3

import subprocess
import configparser
import os

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