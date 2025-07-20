#!/bin/bash
# Download and visualize 6WVR
echo "Downloading 6WVR.cif..."
wget -O 6WVR.cif https://files.rcsb.org/download/6WVR.cif
echo "Starting PyMOL..."
pymol 6WVR.cif -d "@6wvr_visualize.pml"
