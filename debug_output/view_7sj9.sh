#!/bin/bash
# Download and visualize 7SJ9
echo "Downloading 7SJ9.cif..."
wget -O 7SJ9.cif https://files.rcsb.org/download/7SJ9.cif
echo "Starting PyMOL..."
pymol 7SJ9.cif -d "@7sj9_visualize.pml"
