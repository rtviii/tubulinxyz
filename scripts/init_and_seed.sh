#!/bin/sh
# Bootstrap script for fresh deployments.
# Runs as a long-lived service: initializes DB, collects all structures, uploads them.
# Writes live progress to /var/log/tubxz/bootstrap_status.json so the frontend can show it.

STATUS_FILE="/var/log/tubxz/bootstrap_status.json"
mkdir -p /var/log/tubxz

write_status() {
    cat > "$STATUS_FILE" << STATUSEOF
{
  "phase": "$1",
  "message": "$2",
  "collected": $3,
  "uploaded": $4,
  "total_expected": $5,
  "errors": $6,
  "done": $7,
  "timestamp": "$(date -Iseconds)"
}
STATUSEOF
}

write_status "starting" "Initializing database..." 0 0 0 0 false

# Step 1: Initialize DB constraints and phylogeny tree (idempotent)
echo "=== TubXYZ Bootstrap: $(date -Iseconds) ==="
echo "Initializing database constraints and phylogeny..."
python cli.py init-db
write_status "init_db" "Database initialized. Loading literature data..." 0 0 0 0 false

# Step 1.5: Ingest Morisette literature data (mutations + post-translational
# modifications) BEFORE the long structure collection phase. Idempotent (MERGE).
# FATAL on failure: PTMs are core data and a silently-empty Modification table
# is worse than a failed bootstrap (it looks like the app works, but the PTM
# panel is just empty). Source: Morisette et al. 2023
# (doi:10.1371/journal.pone.0295279); CSVs in data/<family>_tubulin/*.csv.
echo "Ingesting Morisette mutation/PTM literature data..."
python -m lib.etl.ingest_morisette --family tubulin_alpha || {
    write_status "error" "Morisette alpha PTM ingest failed" 0 0 0 1 false
    echo "FATAL: Morisette alpha ingest failed -- aborting bootstrap so PTMs cannot be silently missing."
    exit 1
}
python -m lib.etl.ingest_morisette --family tubulin_beta || {
    write_status "error" "Morisette beta PTM ingest failed" 0 0 0 1 false
    echo "FATAL: Morisette beta ingest failed -- aborting bootstrap so PTMs cannot be silently missing."
    exit 1
}

# Post-ingest sanity: ensure each family produced Modification rows AND each
# row carries the species fingerprint we expect (tax_id + OCCURS_IN edge).
# Catches: empty CSV, UTN mapper drops every row, schema regression where
# tax_id stops being written, or PhylogenyNode merge silently fails.
echo "Verifying PTM ingest landed rows + species linkage in Neo4j..."
python -c "
import os, sys
from neo4j import GraphDatabase
drv = GraphDatabase.driver(
    os.environ['NEO4J_URI'],
    auth=(os.environ['NEO4J_USER'], os.environ['NEO4J_PASSWORD']),
)
errors = []
with drv.session(database=os.environ.get('NEO4J_CURRENTDB', 'neo4j')) as s:
    for fam in ('tubulin_alpha', 'tubulin_beta'):
        row = s.run('''
            MATCH (m:Modification {family: \$f})
            RETURN count(m)                                              AS n,
                   sum(CASE WHEN m.tax_id IS NULL THEN 1 ELSE 0 END)      AS no_tax_id,
                   sum(CASE WHEN NOT (m)-[:OCCURS_IN]->(:PhylogenyNode)
                            THEN 1 ELSE 0 END)                            AS no_edge
        ''', f=fam).single()
        n, no_tax_id, no_edge = row['n'], row['no_tax_id'], row['no_edge']
        print(f'  {fam}: {n} mods, missing tax_id: {no_tax_id}, missing OCCURS_IN: {no_edge}')
        if n == 0:
            errors.append(f'{fam}: 0 Modification rows after ingest')
        if no_tax_id:
            errors.append(f'{fam}: {no_tax_id} Modifications without tax_id')
        if no_edge:
            errors.append(f'{fam}: {no_edge} Modifications without OCCURS_IN edge')
drv.close()
if errors:
    for e in errors:
        print(f'FATAL: {e}', file=sys.stderr)
    sys.exit(1)
" || {
    write_status "error" "PTM post-ingest verification failed" 0 0 0 1 false
    exit 1
}

# Step 2: Collect every structure that lacks a COMPLETE local profile.
# Idempotent and resumable:
#   - structures with a complete <id>.json are skipped fast (generate_profile
#     short-circuits on overwrite=False)
#   - structures only partially collected by an interrupted previous run (have
#     a dir + intermediate files but no <id>.json) are regenerated to completion
# This replaces the old "if any .json exists, skip collection" branch, which
# orphaned partially-collected structures when a previous run was interrupted
# (e.g. by `docker compose down` during a redeploy).
echo "Collecting structures missing a complete local profile..."
write_status "collecting" "Collecting structures from PDB..." 0 0 809 0 false
# HANDOFF NOTE (for the session integrating this into the image/scripts):
# Incremental upload -- collect a structure, then immediately ingest it into
# Neo4j, in one sequential pass -- is the intended DEFAULT ingestion method.
# The old "collect everything, then upload everything" two-phase flow should be
# retired across ALL ingestion paths:
#   1. cli.py `collect-missing`: flip `--upload` to default True (or make it
#      unconditional). This script already passes --upload (see below).
#   2. cli.py `weekly-ingest`: it currently runs a collect phase then a separate
#      upload phase -- rework it to interleave the same way (ingest each
#      structure right after collecting it).
#   3. Keep it single-threaded: collect and upload must never run as concurrent
#      Neo4j writers. Do not reintroduce a parallel/background uploader.
# Why: the catalogue fills in as structures are collected, instead of sitting
# empty for the whole multi-hour collection phase.
# --upload makes collect-missing ingest each structure into Neo4j right after
# it is collected, so the catalogue page fills in incrementally during this
# (multi-hour) phase instead of staying empty until step 3. collect-missing is
# single-threaded, so this stays strictly sequential with step 3 below -- the
# two never run a Neo4j writer concurrently.
UPLOADED=0
python cli.py collect-missing --upload 2>&1 | while IFS= read -r line; do
    echo "$line"
    case "$line" in
        *"Uploaded "*" to Neo4j"*) UPLOADED=$((UPLOADED + 1)) ;;
    esac
    COUNT=$(find /mnt/tubetl_data -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
    write_status "collecting" "Collecting structures and filling the catalogue..." "$COUNT" "$UPLOADED" 809 0 false
done

# Step 3: Catch any profiles that are on disk but still missing from Neo4j --
# transient ingest failures during step 2, or profiles left behind by a
# previous interrupted run. Idempotent; usually a fast no-op now that step 2
# uploads as it goes.
echo "Uploading any local profiles still missing from Neo4j..."
write_status "uploading" "Finalizing catalogue..." 0 0 0 0 false
python cli.py upload-missing --workers 4

write_status "done" "Bootstrap complete." 0 0 0 0 true
echo "=== TubXYZ Bootstrap complete: $(date -Iseconds) ==="

# Keep container alive so docker compose doesn't restart it.
# The scheduler handles weekly updates from here.
echo "Bootstrap finished. Sleeping forever (scheduler handles weekly updates)."
tail -f /dev/null
