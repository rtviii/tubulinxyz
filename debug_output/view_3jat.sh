#!/bin/bash
# Download and visualize 3JAT
echo "Downloading 3JAT.cif..."
wget -O 3JAT.cif https://files.rcsb.org/download/3JAT.cif
echo "Starting PyMOL..."
pymol 3JAT.cif -d "@3jat_visualize.pml"
