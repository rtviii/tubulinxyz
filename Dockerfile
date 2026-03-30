# ---- Stage 1: Node.js dependencies (tsx, molstar for CIF extraction) ----
# gl (headless OpenGL) requires python3 + build tools for native compilation
FROM node:20-slim AS node-deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 make g++ pkg-config \
    libxi-dev libgl1-mesa-dev \
    && ln -sf /usr/bin/python3 /usr/bin/python \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY package.json package-lock.json* ./
RUN npm install --production

# ---- Stage 2: Python dependencies (compiled without build-essential in final image) ----
FROM python:3.11-slim AS py-builder
RUN apt-get update && apt-get install -y --no-install-recommends build-essential && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir --prefix=/install --timeout 120 -r requirements.txt

# ---- Stage 3: Final image ----
FROM python:3.11-slim
WORKDIR /app

# Node.js runtime (needed for tsx/molstar CIF extraction in ETL pipeline)
COPY --from=node-deps /usr/local/bin/node /usr/local/bin/node
COPY --from=node-deps /usr/local/lib/node_modules /usr/local/lib/node_modules
RUN ln -s /usr/local/lib/node_modules/npm/bin/npm-cli.js /usr/local/bin/npm
COPY --from=node-deps /app/node_modules ./node_modules

# Python dependencies
COPY --from=py-builder /install /usr/local

# System packages: cron (scheduler), curl (healthchecks)
RUN apt-get update && apt-get install -y --no-install-recommends cron curl && rm -rf /var/lib/apt/lists/*

# Application code
COPY . .
RUN chmod +x /app/muscle3.8.1

# Download NCBI taxonomy database at build time (ete3).
# This eliminates the need to mount ncbi_taxonomy.sqlite from the host.
ENV NCBI_TAXA_SQLITE=/app/.etetoolkit/taxa.sqlite
RUN mkdir -p /app/.etetoolkit && python -c "from ete3 import NCBITaxa; NCBITaxa(dbfile='/app/.etetoolkit/taxa.sqlite')"

# Crontab for the scheduler service (weekly ingestion).
# Only used when the container is started with `cron -f` as entrypoint.
COPY scripts/crontab /etc/cron.d/tubxz-ingest
RUN chmod 0644 /etc/cron.d/tubxz-ingest && crontab /etc/cron.d/tubxz-ingest

EXPOSE 8000
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000", "--root-path", "/api"]
