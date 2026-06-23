#!/usr/bin/env python3
from pathlib import Path

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import configparser
import argparse
import subprocess
import os

from helper import parse_ini, get_db_name


def populated_file(path: str) -> bool:
    """Return True when path exists and has content."""
    return os.path.exists(path) and os.path.getsize(path) > 0


def parse_args(args=None):
    """Parse command-line arguments for generate_chrom_sizes.

    Args:
        args (list|None): Optional argument list for testing.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(dest="species", type=str, help="species production name")
    parser.add_argument(dest="assembly", type=str, help="assembly default")
    parser.add_argument(dest="version", type=int, help="Ensembl release version")
    parser.add_argument(
        "-I",
        "--ini_file",
        dest="ini_file",
        type=str,
        required=False,
        help="full path database configuration file, default - DEFAULT.ini in the same directory.",
    )
    parser.add_argument(
        "--chrom_sizes",
        dest="chrom_sizes",
        type=str,
        required=False,
        help="file with chromomsome sizes, default - <species>_<assembly>.chrom.sizes in the same directory.",
    )
    parser.add_argument(
        "--fasta",
        dest="fasta",
        type=str,
        required=False,
        help="Bgzipped FASTA to use for chrom sizes when no core database is available.",
    )
    parser.add_argument(
        "--force",
        dest="force",
        action="store_true",
        help="forcefully create config even if already exists",
    )

    return parser.parse_args(args)


def generate_chrom_sizes_from_fasta(
    fasta: str, chrom_sizes: str, force: bool = False
) -> None:
    """Generate chromosome sizes from a FASTA .fai index."""
    if populated_file(chrom_sizes) and not force:
        print(f"[INFO] {chrom_sizes} file already exists, skipping ...")
        return

    if not fasta or not os.path.isfile(fasta):
        raise FileNotFoundError(f"Could not generate chrom sizes. FASTA missing - {fasta}")

    fai = fasta + ".fai"
    if not os.path.isfile(fai):
        process = subprocess.run(
            ["samtools", "faidx", fasta],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if process.returncode != 0:
            raise Exception(
                f"Could not index FASTA for chrom sizes - {fasta}\n{process.stderr.decode().strip()}"
            )

    with open(fai, "r") as in_file, open(chrom_sizes, "w") as out_file:
        for line in in_file:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            name, length = parts[0], parts[1]
            out_file.write(f"{name}\t{int(length) + 1}\n")


def generate_chrom_sizes(
    server: dict,
    core_db: str,
    chrom_sizes: str,
    assembly: str = "grch38",
    force: bool = False,
) -> None:
    """Generate a chromosome sizes file from the core database.

    Writes seq_region lengths and synonym lengths with deduplication.

    Args:
        server (dict): Server connection mapping.
        core_db (str): Core database name.
        chrom_sizes (str): Output chrom sizes filename.
        assembly (str): Assembly identifier used to filter coord_system.
        force (bool): If True overwrite existing file.

    Returns:
        None
    """
    if populated_file(chrom_sizes) and not force:
        print(f"[INFO] {chrom_sizes} file already exists, skipping ...")
        return

    query = f"SELECT coord_system_id FROM coord_system WHERE version = '{assembly}';"
    process = subprocess.run(
        [
            "mysql",
            "--host",
            server["host"],
            "--port",
            server["port"],
            "--user",
            server["user"],
            "--database",
            core_db,
            "-N",
            "--execute",
            query,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    coord_ids = (
        "(" + ",".join([id for id in process.stdout.decode().strip().split("\n")]) + ")"
    )

    query = f"SELECT name, length FROM seq_region WHERE coord_system_id IN {coord_ids};"
    process = subprocess.run(
        [
            "mysql",
            "--host",
            server["host"],
            "--port",
            server["port"],
            "--user",
            server["user"],
            "--database",
            core_db,
            "-N",
            "--execute",
            query,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    with open(chrom_sizes, "w") as file:
        file.write(process.stdout.decode())

    query = f"SELECT ss.synonym, s.length FROM seq_region AS s, seq_region_synonym AS ss WHERE s.seq_region_id = ss.seq_region_id AND s.coord_system_id IN {coord_ids};"
    process = subprocess.run(
        [
            "mysql",
            "--host",
            server["host"],
            "--port",
            server["port"],
            "--user",
            server["user"],
            "--database",
            core_db,
            "-N",
            "--execute",
            query,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    with open(chrom_sizes, "a") as file:
        out = process.stdout.decode().strip()
        if out:
            if Path(chrom_sizes).stat().st_size > 0:
                file.write("\n")
            file.write(out)
            file.write("\n")

    # remove duplicates
    with open(chrom_sizes, "r") as file:
        lines = file.readlines()

    lengths = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue
        parts = [col.strip() for col in line.split("\t")]
        if len(parts) != 2:
            sys.exit(f"ERROR: invalid chrom_sizes row: {line!r}")
        name, length = parts
        if name not in lengths or int(lengths[name]) < int(length):
            lengths[name] = length

    with open(chrom_sizes, "w") as file:
        for name in lengths:
            # we will keep length + 1 because bedToBigBed fails if it finds variant at boundary
            length = int(lengths[name]) + 1
            file.write(f"{name}\t{str(length)}\n")


def main(args=None):
    """Main entry point to create chromosome sizes file.

    Parses arguments, determines the core DB and invokes generate_chrom_sizes.

    Args:
        args (list|None): Optional argument list for testing.

    Returns:
        None
    """
    args = parse_args(args)

    species = args.species
    assembly = args.assembly
    chrom_sizes = args.chrom_sizes or f"{species}_{assembly}.chrom.sizes"
    ini_file = args.ini_file or "DEFAULT.ini"
    core_server = parse_ini(ini_file, "core")
    core_db = get_db_name(core_server, args.version, species, type="core")
    if core_db == "" or core_db is None:
        print(
            f"[WARNING] Could not resolve core database for {species} release {args.version}; generating chrom sizes from FASTA"
        )
        generate_chrom_sizes_from_fasta(args.fasta, chrom_sizes, args.force)
        return

    generate_chrom_sizes(core_server, core_db, chrom_sizes, assembly, args.force)
    if not populated_file(chrom_sizes):
        print(
            f"[WARNING] Core database {core_db} did not produce chrom sizes for {species} {assembly}; generating chrom sizes from FASTA"
        )
        generate_chrom_sizes_from_fasta(args.fasta, chrom_sizes, True)


if __name__ == "__main__":
    sys.exit(main())
