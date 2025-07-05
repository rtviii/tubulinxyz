# PyMOL script for 6E7B N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/6E7B.cif)
load 6E7B.cif
hide everything
show cartoon
color grey, all

# Color N-terminus residues in green
select nterm_B, chain B and resi 1
color green, nterm_B
show spheres, nterm_B
select nterm_A, chain A and resi 1
color green, nterm_A
show spheres, nterm_A

# Show connections as distances
distance conn_B_A, nterm_B, nterm_A

# Color connected pairs in red
color red, chain B or chain A

# Summary: 1 connections from 2 chains
# Success rate: 50.0%

zoom all
