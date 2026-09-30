# Phenotypes pipeline

Pipeline to download, parse, clean and load phenotype data into PostgreSQL.
This pipeline uses Nextflow 26.04.0.

```bash
module load nextflow/26.04.0
```

## Project setup

Run commands in this README from the project directory:

```bash
cd ensembl-variation-pipelines/nextflow/phenotypes
```

The Python project requires Python 3.14 and uses `uv`. Create or update the virtual environment with:

```bash
uv sync
```

## Database connection

Copy the example connection file:

```bash
cp .env.example .env
```

Set these values in `.env`:

```.env
PGDATABASE=database_name
PGUSER=database_user
PGHOST=database_host
PGPORT=database_port
PGSCHEMA=phenotypes
```

Load `.env`:

```bash
set -a; source .env; set +a
```

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
./bin/apply_schema.sh
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

## ClinVar

By default, the pipeline downloads the current [ClinVar RCV XML release](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_00-latest.xml.gz).
Failed download attempts are retried three times.

To use a manually downloaded ClinVar XML file instead:

```bash
nextflow run main.nf -profile standard --input /path/to/ClinVarRCVRelease_00-latest.xml.gz
```

Using `--input` skips the download.

The Python parser is currently run separately from Nextflow. It retains supported
source associations even when they have no Ensembl match, including single-measure
structural variants. A source placement on the requested assembly is required:
records with no such placement are skipped as `missing_requested_assembly_location`,
even if they have an rsID. Imprecise placements on that assembly remain supported.
Structural status produces a warning, not a rejection. Haplotypes and other
unsupported compound representations are still rejected explicitly.

Source coordinates are stored once in
`source_report.reported_variant.locations`, including uncertain bounds and the
source VCF representation when supplied. Only the requested assembly is retained;
unknown coordinates stay null. Source variant type and copy-number/cytogenetic
attributes are also preserved. The parser does not emit empty resolution
placeholders; the resolver will add its results, provide optional website links and flag
ambiguous/conflicting matches for review; absence of a match will not exclude a
supported association from the explorer dataset.

`clinvar_output.jsonl` contains data for loading and variant lookup, without duplicated
accessions, trait IDs/types, record statuses or sample origins. The RCV accession
is stored in `source_report.source_accession`; use it to find the original XML
record when debugging. The release date is stored once in `clinvar_summary.json`,
which should accompany the records through resolution and loading for `ImportRun`.
For easy checking, the parser writes two tab-separated files beside the summary:
`warnings.txt` with columns
`RCV`, `VCV`, `warning`, `details`, and `rejections.txt` with columns `RCV`, `VCV`,
`rejection_reason`, `details`. Each warning occurrence or rejected RCV gets a row.
Rejection details contain the source value, including requested/available assemblies;
structured values are JSON within the column. Details are blank when unavailable,
as are missing accessions. Headers are always written, even for empty reports.

`clinvar_summary.json` records input RCVs (`n_records_seen`), accepted
RCVs and rejected RCVs, with percentages of the input total rounded to two decimal
places. `rejections_by_reason` gives the rejected count for each reason.
`warnings` gives the occurrence count for each warning type; one RCV can produce
multiple warnings, so these are not counts of distinct warned RCVs.
`n_records_emitted` counts JSONL rows separately, because one RCV can produce
several rows. `n_unique_phenotypes` counts distinct output phenotype names, ignoring
case and repeated whitespace; this is not ontology-based grouping. Metadata includes
`source`, `release_date` and `requested_assembly`. Empty input gives zero percentages.
Fields are written in this order: metadata, input count, accepted/output counts,
rejected counts and reasons, then warnings. No by-level assessment counts are included.
Summary count keys have an `n_` prefix and are source-neutral: `n_records_seen`,
`n_records_accepted` and `n_records_rejected` count input source records (RCVs for ClinVar), while
`n_records_emitted` counts output JSONL rows.
