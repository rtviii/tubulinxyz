#!/bin/bash
# Download and visualize 4O2B
echo "Downloading 4O2B.cif..."
wget -O 4O2B.cif https://files.rcsb.org/download/4O2B.cif
echo "Starting PyMOL..."
pymol 4O2B.cif -d "@4o2b_visualize.pml"
