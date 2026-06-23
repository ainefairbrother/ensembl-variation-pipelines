#!/usr/bin/env python3

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
import argparse
import subprocess
import os
import re

from helper import *

GFF_FASTA_BASE_DIR = (
    "/hps/nobackup/flicek/ensembl/production/ensembl_dumps/ftp_mvp/organisms"
)
GFF_FASTA_BASE_URL = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"
GFF_FILE_NAME = "sorted_genes.gff3.gz"
REMOTE_GFF_FILE_NAME = "genes.gff3.gz"


def format_geneset_date(value):
    if not value:
        return None
    return re.sub("[-\\s]", "_", value)


def get_genome_genebuild_date(server, metadata_db, genome_uuid):
    query = f"SELECT genebuild_date FROM genome WHERE genome_uuid = '{genome_uuid}';"
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
            metadata_db,
            "-N",
            "--execute",
            query,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process.returncode != 0:
        print(f"[WARNING] Failed to retrieve genome.genebuild_date for genome - {genome_uuid}")
        print(f"\tError - {process.stderr.decode().strip()}")
        return None
    return process.stdout.decode().strip() or None


def parse_args(args=None):
    """Parse command-line arguments for processing a GFF.

    Args:
        args (list|None): Optional argument list for testing.

    Returns:
        argparse.Namespace: Parsed arguments including genome_uuid, release_id, out_dir,
            ini_file, gff_dir and force flag.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="genome_uuid", type=str, help="Genome uuid")
    parser.add_argument(
        dest="release_id", type=str, help="Ensembl release id from metadata database"
    )
    parser.add_argument(
        "--out_dir",
        dest="out_dir",
        type=str,
        help="Out directory where processed GFF file will be created",
    )
    parser.add_argument(
        "-I",
        "--ini_file",
        dest="ini_file",
        type=str,
        required=False,
        help="Full path database configuration file, default - DEFAULT.ini in the same directory.",
    )
    parser.add_argument(
        "--gff_dir", dest="gff_dir", type=str, required=False, help="GFF directory"
    )
    parser.add_argument(
        "--ftp_cache_dir",
        dest="ftp_cache_dir",
        type=str,
        required=False,
        help="Writable local cache root mirroring pub/ensemblorganisms for GFF downloads",
    )
    parser.add_argument(
        "--download_retries",
        dest="download_retries",
        type=int,
        default=3,
        help="Number of download attempts for remote FTP fallback",
    )
    parser.add_argument("--force", dest="force", action="store_true")

    return parser.parse_args(args)


def index_gff(bgzipped_gff: str, force: str = False) -> None:
    """Create a CSI index for a bgzipped GFF file using tabix.

    Args:
        bgzipped_gff (str): Path to the bgzipped GFF file.
        force (bool): If False and a .csi exists, skip indexing.

    Raises:
        FileNotFoundError: If the GFF file does not exist.
        SystemExit: Exits with error code 1 if tabix indexing fails.
    """
    if not os.path.isfile(bgzipped_gff):
        raise FileNotFoundError(
            f"Could not run tabix index. File does not exist - {bgzipped_gff}"
        )

    csi = bgzipped_gff + ".csi"

    if os.path.isfile(csi) and not force:
        print(f"[INFO] {csi} file exist. Skipping ...")
        return

    process = subprocess.run(
        ["tabix", "-f", "-C", bgzipped_gff],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process.returncode != 0:
        print(
            f"[ERROR] Cannot index - {bgzipped_gff}\n{process.stderr.decode()}\nExiting ..."
        )
        exit(1)


def sort_gff(file: str, sorted_file: str = None) -> str:
    """Sort a GFF file by seqname and start/end positions while preserving header lines.

    Args:
        file (str): Path to the input GFF file.
        sorted_file (str|None): Optional output path for the sorted file. If omitted,
            a file named "sorted_<basename>" is created in the same directory.

    Returns:
        str: Path to the sorted file.

    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"Could not sort. File does not exist - {file}")

    sorted_file = sorted_file or os.path.join(
        os.path.dirname(file), "sorted_" + os.path.basename(file)
    )

    os.system(
        f"(grep '^#' {file} & grep -v '^#' {file} | sort -k1,1 -k4,4n -k5,5n -t$'\\t') > {sorted_file}"
    )

    return sorted_file


def main(args=None):
    """Main entry point for processing and indexing a GFF for a genome.

    Uses metadata to locate the appropriate GFF source when necessary, copies, sorts,
    bgzips and indexes the GFF file.

    Args:
        args (list|None): Optional argument list for testing; if None uses sys.argv.

    Returns:
        None
    """
    args = parse_args(args)

    genome_uuid = args.genome_uuid
    release_id = args.release_id
    out_dir = args.out_dir or os.getcwd()
    ini_file = args.ini_file or "DEFAULT.ini"
    gff_dir = args.gff_dir or os.getcwd()
    metadb_server = parse_ini(ini_file, "metadata")

    source_gff = os.path.join(gff_dir, GFF_FILE_NAME)

    if (
        not os.path.isfile(source_gff)
        or not os.path.isfile(source_gff + ".csi")
        or args.force
    ):
        scientific_name = get_scientific_name(
            metadb_server, "ensembl_genome_metadata", genome_uuid
        ).replace(" ", "_")
        if scientific_name == "" or scientific_name is None:
            raise Exception(
                f"[ERROR] Could not retrieve scientific name for genome uuid - {genome_uuid}"
            )
        scientific_name = re.sub("[^a-zA-Z0-9]+", " ", scientific_name)
        scientific_name = re.sub(" +", "_", scientific_name)
        scientific_name = re.sub("^_+|_+$", "", scientific_name)
        assembly_accession = get_assembly_accession(
            metadb_server, "ensembl_genome_metadata", genome_uuid
        )
        if assembly_accession == "" or assembly_accession is None:
            raise Exception(
                f"[ERROR] Could not retrieve assembly accession for genome uuid - {genome_uuid}"
            )
        annotation_source = get_dataset_attribute_value(
            metadb_server,
            "ensembl_genome_metadata",
            genome_uuid,
            release_id,
            "genebuild.annotation_source",
        )
        if annotation_source == "" or annotation_source is None:
            annotation_source = "ensembl"
            print(
                f"[WARNING] Could not retrieve genebuild annotation source for genome uuid - {genome_uuid} and release id - {release_id}; using '{annotation_source}'"
            )
        else:
            annotation_source = annotation_source.lower()
        last_geneset_update = get_dataset_attribute_value(
            metadb_server,
            "ensembl_genome_metadata",
            genome_uuid,
            release_id,
            "genebuild.last_geneset_update",
        )
        if last_geneset_update:
            last_geneset_update = format_geneset_date(last_geneset_update)
        else:
            genome_genebuild_date = get_genome_genebuild_date(
                metadb_server, "ensembl_genome_metadata", genome_uuid
            )
            if genome_genebuild_date:
                last_geneset_update = format_geneset_date(genome_genebuild_date)
                print(
                    f"[WARNING] Could not retrieve last genebuild update date for genome uuid - {genome_uuid} and release id - {release_id}; using genome.genebuild_date {last_geneset_update}"
                )
            else:
                last_geneset_update = None
                print(
                    f"[WARNING] Could not retrieve last genebuild update date for genome uuid - {genome_uuid} and release id - {release_id}; using latest available geneset directory"
                )

        source_gff = find_local_mvp_file(
            GFF_FASTA_BASE_DIR,
            scientific_name,
            assembly_accession,
            REMOTE_GFF_FILE_NAME,
            annotation_source,
            last_geneset_update,
        )

        if not source_gff or not os.path.isfile(source_gff):
            source_gff = find_local_mvp_file(
                args.ftp_cache_dir,
                scientific_name,
                assembly_accession,
                REMOTE_GFF_FILE_NAME,
                annotation_source,
                last_geneset_update,
            )

        if source_gff:
            print(f"[INFO] Using GFF - {source_gff}")
        elif args.ftp_cache_dir:
            cache_root = os.path.join(
                args.ftp_cache_dir,
                scientific_name,
                assembly_accession,
                annotation_source,
                "geneset",
            )
            os.makedirs(cache_root, exist_ok=True)
            downloaded = False
            for url in remote_mvp_candidate_urls(
                GFF_FASTA_BASE_URL,
                scientific_name,
                assembly_accession,
                REMOTE_GFF_FILE_NAME,
                annotation_source,
                last_geneset_update,
            ):
                date_dir = url.rstrip("/").split("/")[-2]
                cache_dir = (
                    os.path.join(cache_root, date_dir)
                    if date_dir != "geneset"
                    else cache_root
                )
                os.makedirs(cache_dir, exist_ok=True)
                source_gff = os.path.join(cache_dir, REMOTE_GFF_FILE_NAME)
                returncode = download_file_with_retries(
                    source_gff, url, args.download_retries
                )
                if returncode == 0:
                    downloaded = True
                    break
            if not downloaded:
                raise FileNotFoundError(
                    f"Could not find GFF locally or download from {GFF_FASTA_BASE_URL} for {assembly_accession}"
                )
        else:
            raise FileNotFoundError(
                f"Could not find - {source_gff}; provide --ftp_cache_dir to download from remote FTP"
            )

        compressed_gff = os.path.join(out_dir, REMOTE_GFF_FILE_NAME)
        returncode = copyto(source_gff, compressed_gff)
        if returncode != 0:
            raise Exception(
                f"Failed to copy.\n\tSource - {source_gff}\n\tTarget - {compressed_gff}"
            )

        unzipped_gff = ungzip_file(compressed_gff)
        sorted_gff = sort_gff(unzipped_gff)
        bgzipped_gff = bgzip_file(sorted_gff)
        index_gff(bgzipped_gff, force=args.force)


if __name__ == "__main__":
    sys.exit(main())
