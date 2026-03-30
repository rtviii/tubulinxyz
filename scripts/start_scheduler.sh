#!/bin/sh
# Entrypoint for the scheduler container.
# Dumps the container environment to a file that cron jobs can source,
# then starts cron in the foreground.

# Dump all env vars so cron jobs can source them
env | grep -E '^(NEO4J_|TUBETL_|TUBXZ_|NCBI_|PATH=|HOME=|PYTHONPATH=)' > /app/.env.cron

mkdir -p /var/log/tubxz

echo "Scheduler started at $(date -Iseconds). Cron jobs will run weekly (Sunday 3am UTC)."
exec cron -f
