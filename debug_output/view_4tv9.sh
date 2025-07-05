#!/bin/bash
# Download and visualize 4TV9
echo "Downloading 4TV9.cif..."
wget -O 4TV9.cif https://files.rcsb.org/download/4TV9.cif
echo "Starting PyMOL..."
pymol 4TV9.cif -d "@4tv9_visualize.pml"
