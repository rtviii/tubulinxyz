#!/bin/bash
# Batch script to download PDB and visualize N-terminus analysis

PDB_ID="6O2T"
DEBUG_DIR="debug_output"

echo "🔽 Downloading ${PDB_ID}.cif..."
wget -O "${DEBUG_DIR}/${PDB_ID}.cif" "https://files.rcsb.org/download/${PDB_ID}.cif"

echo "🎨 Starting PyMOL visualization..."
cd "${DEBUG_DIR}"
pymol "${PDB_ID}.cif" -d "@6o2t_visualize.pml"

echo "✅ Done! PyMOL should now be showing the N-terminus analysis."
