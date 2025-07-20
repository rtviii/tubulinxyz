#!/bin/bash
# Download and visualize 7SJ8
echo "Downloading 7SJ8.cif..."
wget -O 7SJ8.cif https://files.rcsb.org/download/7SJ8.cif
echo "Starting PyMOL..."
pymol 7SJ8.cif -d "@7sj8_visualize.pml"
