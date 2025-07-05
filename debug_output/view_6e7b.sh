#!/bin/bash
# Download and visualize 6E7B
echo "Downloading 6E7B.cif..."
wget -O 6E7B.cif https://files.rcsb.org/download/6E7B.cif
echo "Starting PyMOL..."
pymol 6E7B.cif -d "@6e7b_visualize.pml"
