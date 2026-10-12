# Phenotypes pipeline

Pipeline to download, parse, clean and load phenotype data into PostgreSQL.
This pipeline uses Nextflow 26.04.0.

```bash
module load nextflow/26.04.0
```

## Project setup

Create the Python environment in the project directory:

```bash
cd ensembl-variation-pipelines/nextflow/phenotypes
PIPELINE_DIR="$PWD"
```

The Python project requires Python 3.14 and uses `uv`. Create or update the virtual environment with:

```bash
uv sync
```

## Database connection

Copy `.env` to the run directory:

```bash
cd /path/to/run-directory
cp -n "$PIPELINE_DIR/.env.example" .env
chmod 600 .env
```

Set the postgres db connection details in `.env`:

```.env
[postgresql]
host=database_host
port=database_port
user=database_user
name=database_name
schema=phenotypes
password=
```

# First time only set-up **(first time only)**

## PostgreSQL setup **(first time only)**

The schema requires PostgreSQL 10 or later and was created using PostgreSQL 16.2. Load the `postgresql/16` module:

```bash
module load postgresql/16
```

Create the schema:

```sql
CREATE SCHEMA phenotypes AUTHORIZATION ensevp;
```

## Create the tables **(first time only)**

Apply `sql/schema.sql` to the configured empty schema:

```bash
"$PIPELINE_DIR/bin/apply_schema.sh"
```

Check tables got generated correctly:

```bash
psql
```

At psql prompt: 

```SQL
-- list schemas
\dn

-- select phenotype schema
SET search_path TO phenotypes;

-- list tables
\dt
```

# Data sources

## ClinVar

By default, the pipeline downloads the current [ClinVar RCV XML release](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_00-latest.xml.gz)
and [VCV XML release](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarVCVRelease_00-latest.xml.gz).
Failed download attempts are retried three times.

To use manually downloaded files, set `sources.clinvar.inputs.rcv` and/or `sources.clinvar.inputs.vcv` in the source configuration (`assets/source_config.json`). Inputs, species and assembly also come from `assets/source_config.json`. At present, these are hard-coded and can be edited by the user running the pipeline. 

By default, the source config file below is supplied to the pipeline. Here, `rcv` and `vcv` are `null`, and so the pipeline will download the files for you. You can supply an alternative `source_config.json` file using `nextflow run main.nf" --source_config /path/to/source_config.json`

```json
{
  "sources": {
    "clinvar": {
      "inputs": {
        "rcv": null,
        "vcv": null
      },
      "contexts": [
        {
          "species_taxonomy_id": 9606,
          "species_production_name": "homo_sapiens",
          "assembly": "GRCh38",
          "assembly_uuid": "fd7fea38-981a-4d73-a879-6f9daef86f08"
        }
      ]
    }
  }
}
```