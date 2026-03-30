#!/bin/sh
set -e

echo "=== TubXYZ Init: $(date -Iseconds) ==="

# Step 1: Initialize DB constraints and phylogeny tree (idempotent)
echo "Initializing database constraints and phylogeny..."
python cli.py init-db

# Step 2: If TUBETL_DATA is empty, do a full collection + upload.
# This handles the first-deploy case where no data exists yet.
# On subsequent deploys (e.g. after docker compose down/up), this is a no-op
# because the data volume already has profiles.
PROFILE_COUNT=$(find /mnt/tubetl_data -maxdepth 2 -name "*.json" -type f 2>/dev/null | head -5 | wc -l)

if [ "$PROFILE_COUNT" -eq 0 ]; then
    echo "TUBETL_DATA is empty -- running full collection (this will take 2-3 hours)..."
    python cli.py collect-all --workers 4
    echo "Uploading all collected profiles to Neo4j..."
    python cli.py upload-all --workers 4
else
    echo "TUBETL_DATA already has data ($PROFILE_COUNT+ profiles found). Skipping bootstrap collection."
    echo "Uploading any missing profiles to Neo4j..."
    python cli.py upload-missing --workers 4
fi

echo "=== TubXYZ Init complete: $(date -Iseconds) ==="
