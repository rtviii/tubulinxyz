#!/bin/sh
# Nuke the local dev DB and rebuild it from scratch, mirroring what the deploy's
# bootstrap (scripts/init_and_seed.sh) does: fresh Community container + empty
# data volume, schema + phylogeny, Morisette mutations/PTMs, then every structure
# profile already on disk under $TUBETL_DATA (no PDB collection).
#
#   sh scripts/dev_reset.sh
set -e
cd "$(dirname "$0")/.."

. ./.env

PY=./venv/bin/python
[ -x "$PY" ] || PY=python

echo "== Removing dev container and data volume"
docker rm -fv tubxz-dev-neo4j >/dev/null 2>&1 || true
docker volume rm tubxz-dev-neo4j-data >/dev/null 2>&1 || true

echo "== Starting fresh Neo4j"
sh scripts/dev_neo4j.sh

echo "== Schema + phylogeny"
"$PY" cli.py init-db

echo "== Morisette mutations/PTMs"
"$PY" -m lib.etl.ingest_morisette --family tubulin_alpha
"$PY" -m lib.etl.ingest_morisette --family tubulin_beta

echo "== Uploading local profiles from $TUBETL_DATA"
"$PY" cli.py upload-missing --workers 4

echo "== Dev DB rebuilt."
