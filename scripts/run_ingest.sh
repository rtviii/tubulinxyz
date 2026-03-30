#!/bin/sh
# Wrapper for cron: loads container environment variables (cron doesn't inherit them)
# then runs the weekly-ingest CLI command.

# Import environment from the file dumped at container start
if [ -f /app/.env.cron ]; then
    . /app/.env.cron
fi

mkdir -p /var/log/tubxz

echo "=== Cron ingest started: $(date -Iseconds) ==="
cd /app && /usr/local/bin/python cli.py weekly-ingest --workers 4
echo "=== Cron ingest finished: $(date -Iseconds) ==="
