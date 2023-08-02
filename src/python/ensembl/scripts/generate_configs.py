import configparser
import argparse
import subprocess
import os

# define arguemnts
parser = argparse.ArgumentParser()
parser.add_argument('-I', '--ini_file', dest="ini_file", type=str, required = False, \
                    help="""
                    Provide ini file location from where database configuration will be read.
                    By default it tries to read DEFAULT.ini in the same directory.
                    """)
parser.add_argument('--genome', dest="genome", type=str, required = True, help="Genome string with format <species>_<assembly>")
parser.add_argument('--version', dest="version", type=str, required = True, help="Ensembl release version to use")
parser.add_argument('--force', dest="force", action="store_true", help="Ensembl release version to use")
opts = parser.parse_args()

# define global variables
ini_file = opts.ini_file or "DEFAULT.ini"
genome = opts.genome
version = opts.version
force = opts.force
assembly = genome.split("_")[-1]
species = genome.replace(f"_{assembly}", "")
config_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../../../nextflow/nf_config"

print("current wd: " + config_dir)

# parse ini config file
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

# get the database name 
query = f"SHOW DATABASES LIKE \"{species}%core%{version}%\";"
process = subprocess.run(["mysql",
        "--host", host,
        "--port", port,
        "--user", user,
        "-N",
        "--execute", query
    ],
    stdout = subprocess.PIPE,
    stderr = subprocess.PIPE
)
dbname = process.stdout.decode().strip()

def generate_synonyms():
    synonym_file = f"{config_dir}/synonyms/{genome}.txt"
    
    if os.path.exists(synonym_file) and not force:
        print(f"[INFO] {synonym_file} file already exists, skipping ...")
        return
        
    query = f"SELECT ss.synonym, sr.name FROM seq_region AS sr, seq_region_synonym AS ss WHERE sr.seq_region_id = ss.seq_region_id;"
    process = subprocess.run(["mysql",
            "--host", host,
            "--port", port,
            "--user", user,
            "--database", dbname,
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
    
def generate_chrom_sizes():
    chrom_sizes_file = f"{config_dir}/chrom_sizes/{genome}.chrom.sizes"
    
    if os.path.exists(chrom_sizes_file) and not force:
        print(f"[INFO] {chrom_sizes_file} file already exists, skipping ...")
        return
        
    query = f"SELECT coord_system_id FROM coord_system WHERE version = '{assembly}';"
    process = subprocess.run(["mysql",
            "--host", host,
            "--port", port,
            "--user", user,
            "--database", dbname,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    coord_ids = tuple([int(id) for id in process.stdout.decode().strip().split("\n")])
    
    query = f"SELECT name, length FROM seq_region WHERE coord_system_id IN {coord_ids};"
    process = subprocess.run(["mysql",
            "--host", host,
            "--port", port,
            "--user", user,
            "--database", dbname,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    with open(chrom_sizes_file, "w") as file: 
        file.write(process.stdout.decode())
        
    query = f"SELECT ss.synonym, s.length FROM seq_region AS s, seq_region_synonym AS ss WHERE s.seq_region_id = ss.seq_region_id AND s.coord_system_id IN {coord_ids};"
    process = subprocess.run(["mysql",
            "--host", host,
            "--port", port,
            "--user", user,
            "--database", dbname,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    with open(chrom_sizes_file, "a") as file: 
        file.write(process.stdout.decode().strip())
    
    # remove duplicates    
    with open(chrom_sizes_file, "r") as file: 
        lines = file.readlines()
        
    lengths = {}
    for line in lines:
        name, length = [col.strip() for col in line.split("\t")]  
        if name not in lengths or int(lengths[name]) < int(length): 
            lengths[name] = length
            
    with open(chrom_sizes_file, "w") as file:
        for name in lengths:
            # we will keep length + 1 because bedToBigBed fails if it finds variant at boundary
            length = int(lengths[name]) + 1
            file.write(f"{name}\t{str(length)}\n")
            
generate_synonyms()
generate_chrom_sizes()