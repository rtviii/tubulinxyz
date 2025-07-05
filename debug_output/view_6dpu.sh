#!/bin/bash
# Download and visualize 6DPU
echo "Downloading 6DPU.cif..."
wget -O 6DPU.cif https://files.rcsb.org/download/6DPU.cif
echo "Starting PyMOL..."
pymol 6DPU.cif -d "@6DPU_visualize.pml"
