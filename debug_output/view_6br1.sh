#!/bin/bash
# Download and visualize 6BR1
echo "Downloading 6BR1.cif..."
wget -O 6BR1.cif https://files.rcsb.org/download/6BR1.cif
echo "Starting PyMOL..."
pymol 6BR1.cif -d "@6br1_visualize.pml"
