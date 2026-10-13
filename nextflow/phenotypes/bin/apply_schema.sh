#!/bin/bash

set -euo pipefail

PIPELINE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${1:-$PWD/.env}"
PYTHON_BIN="${PIPELINE_DIR}/.venv/bin/python"
SCHEMA_FILE="${PIPELINE_DIR}/sql/schema.sql"

if [[ ! -x "$PYTHON_BIN" ]]; then
    echo "Python environment not found. Run uv sync in $PIPELINE_DIR." >&2
    exit 1
fi

connection_env=$(PYTHONPATH="${PIPELINE_DIR}/src" "$PYTHON_BIN" -m phenotypes.config \
    --env-file "$ENV_FILE" --postgresql-env)
eval "$connection_env"

psql -X -v ON_ERROR_STOP=1 -v target_schema="$PGSCHEMA" \
    -v schema_file="$SCHEMA_FILE" <<'SQL'
SET search_path TO :"target_schema";
\i :schema_file
SQL
