#!/bin/sh
# Local dev Neo4j: the same Community image the deploy runs (tubxz_deploy
# docker-compose.yml), so nothing Enterprise-only can creep into the code.
#
# Reads NEO4J_PASSWORD from ./.env. Bolt and the browser UI are bound to
# localhost only. Returns once bolt accepts connections.
#
# Data lives in a Docker named volume, not a host bind mount: on Docker Desktop
# (macOS) bind-mounted Neo4j data dirs are flaky -- chown warnings, stale
# handles after delete+recreate, containers dying at startup.
#
#   sh scripts/dev_neo4j.sh          # start (creates the container on first run)
#   docker stop tubxz-dev-neo4j      # stop; data persists in the volume
#   sh scripts/dev_reset.sh          # nuke and rebuild the whole dev DB
set -e
cd "$(dirname "$0")/.."

. ./.env

NAME=tubxz-dev-neo4j
VOLUME=tubxz-dev-neo4j-data
IMAGE=neo4j:5.26-community

if docker container inspect "$NAME" >/dev/null 2>&1; then
    docker start "$NAME" >/dev/null
else
    docker run -d --name "$NAME" \
        -p 127.0.0.1:7687:7687 -p 127.0.0.1:7474:7474 \
        -e NEO4J_AUTH="neo4j/$NEO4J_PASSWORD" \
        -e NEO4J_server_memory_heap_max__size=2G \
        -e NEO4J_server_memory_pagecache_size=1G \
        -v "$VOLUME":/data \
        "$IMAGE" >/dev/null
fi

printf "Waiting for Neo4j"
for _ in $(seq 1 60); do
    if docker exec "$NAME" cypher-shell -u neo4j -p "$NEO4J_PASSWORD" "RETURN 1" >/dev/null 2>&1; then
        echo " ready on bolt://localhost:7687"
        exit 0
    fi
    if [ "$(docker inspect -f '{{.State.Running}}' "$NAME")" != "true" ]; then
        echo " container exited:"
        docker logs --tail 30 "$NAME"
        exit 1
    fi
    printf "."
    sleep 3
done
echo " timed out after 180s"
docker logs --tail 30 "$NAME"
exit 1
