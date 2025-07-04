#!/bin/bash
# Download and visualize 6O2T
echo "Downloading 6O2T.cif..."
wget -O 6O2T.cif https://files.rcsb.org/download/6O2T.cif
echo "Starting PyMOL..."
pymol 6O2T.cif -d "@6o2t_visualize.pml"
