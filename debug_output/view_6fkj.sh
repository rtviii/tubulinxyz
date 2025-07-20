#!/bin/bash
# Download and visualize 6FKJ
echo "Downloading 6FKJ.cif..."
wget -O 6FKJ.cif https://files.rcsb.org/download/6FKJ.cif
echo "Starting PyMOL..."
pymol 6FKJ.cif -d "@6fkj_visualize.pml"
