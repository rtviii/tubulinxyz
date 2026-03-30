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

# Step 2: Figure out what needs to happen
PROFILE_COUNT=$(find /mnt/tubetl_data -maxdepth 2 -name "*.json" -type f 2>/dev/null | head -5 | wc -l)

if [ "$PROFILE_COUNT" -eq 0 ]; then
    # Fresh deploy: collect everything
    echo "TUBETL_DATA is empty -- running full collection..."
    write_status "collecting" "Downloading structures from PDB (this takes 2-3 hours)..." 0 0 809 0 false
    python cli.py collect-missing 2>&1 | while IFS= read -r line; do
        echo "$line"
        # Update count periodically
        COUNT=$(find /mnt/tubetl_data -maxdepth 2 -name "*.json" -type f 2>/dev/null | wc -l)
        write_status "collecting" "Downloading structures from PDB..." "$COUNT" 0 809 0 false
    done

    echo "Uploading all collected profiles to Neo4j..."
    write_status "uploading" "Uploading structures to database..." 0 0 0 0 false
    python cli.py upload-all --workers 4
else
    echo "TUBETL_DATA has data ($PROFILE_COUNT+ profiles). Uploading any missing..."
    write_status "uploading" "Syncing structures to database..." 0 0 0 0 false
    python cli.py upload-missing --workers 4
fi

write_status "done" "Bootstrap complete." 0 0 0 0 true
echo "=== TubXYZ Bootstrap complete: $(date -Iseconds) ==="

# Keep container alive so docker compose doesn't restart it.
# The scheduler handles weekly updates from here.
echo "Bootstrap finished. Sleeping forever (scheduler handles weekly updates)."
tail -f /dev/null
