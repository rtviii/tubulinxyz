#!/bin/bash
# Download and visualize 8VRK
echo "Downloading 8VRK.cif..."
wget -O 8VRK.cif https://files.rcsb.org/download/8VRK.cif
echo "Starting PyMOL..."
pymol 8VRK.cif -d "@8vrk_visualize.pml"
