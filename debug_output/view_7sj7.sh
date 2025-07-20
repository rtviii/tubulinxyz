#!/bin/bash
# Download and visualize 7SJ7
echo "Downloading 7SJ7.cif..."
wget -O 7SJ7.cif https://files.rcsb.org/download/7SJ7.cif
echo "Starting PyMOL..."
pymol 7SJ7.cif -d "@7sj7_visualize.pml"
