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
write_status "init_db" "Database initialized. Checking for structures..." 0 0 0 0 false

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
python cli.py collect-missing 2>&1 | while IFS= read -r line; do
    echo "$line"
    COUNT=$(find /mnt/tubetl_data -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
    write_status "collecting" "Collecting structures from PDB..." "$COUNT" 0 809 0 false
done

# Step 3: Upload every local profile not yet in the database.
# Idempotent: uploads (local complete profiles - structures already in Neo4j).
echo "Uploading local profiles missing from Neo4j..."
write_status "uploading" "Uploading structures to database..." 0 0 0 0 false
python cli.py upload-missing --workers 4

write_status "done" "Bootstrap complete." 0 0 0 0 true
echo "=== TubXYZ Bootstrap complete: $(date -Iseconds) ==="

# Keep container alive so docker compose doesn't restart it.
# The scheduler handles weekly updates from here.
echo "Bootstrap finished. Sleeping forever (scheduler handles weekly updates)."
tail -f /dev/null
