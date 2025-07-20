#!/bin/bash
# Download and visualize 1SA0
echo "Downloading 1SA0.cif..."
wget -O 1SA0.cif https://files.rcsb.org/download/1SA0.cif
echo "Starting PyMOL..."
pymol 1SA0.cif -d "@1sa0_visualize.pml"
