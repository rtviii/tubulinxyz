#!/bin/bash
# Download and visualize 6O2R
echo "Downloading 6O2R.cif..."
wget -O 6O2R.cif https://files.rcsb.org/download/6O2R.cif
echo "Starting PyMOL..."
pymol 6O2R.cif -d "@6o2r_visualize.pml"
