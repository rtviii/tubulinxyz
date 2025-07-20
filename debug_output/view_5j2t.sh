#!/bin/bash
# Download and visualize 5J2T
echo "Downloading 5J2T.cif..."
wget -O 5J2T.cif https://files.rcsb.org/download/5J2T.cif
echo "Starting PyMOL..."
pymol 5J2T.cif -d "@5j2t_visualize.pml"
