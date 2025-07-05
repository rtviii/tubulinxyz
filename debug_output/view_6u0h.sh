#!/bin/bash
# Download and visualize 6U0H
echo "Downloading 6U0H.cif..."
wget -O 6U0H.cif https://files.rcsb.org/download/6U0H.cif
echo "Starting PyMOL..."
pymol 6U0H.cif -d "@6u0h_visualize.pml"
