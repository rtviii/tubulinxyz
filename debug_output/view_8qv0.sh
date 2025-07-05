#!/bin/bash
# Download and visualize 8QV0
echo "Downloading 8QV0.cif..."
wget -O 8QV0.cif https://files.rcsb.org/download/8QV0.cif
echo "Starting PyMOL..."
pymol 8QV0.cif -d "@8qv0_visualize.pml"
