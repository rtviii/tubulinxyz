# PyMOL script for 6WVR N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/6WVR.cif)
load 6WVR.cif
hide everything
show cartoon
color grey, all

# Color N-terminus residues in green
select nterm_A, chain A and resi 1
color green, nterm_A
show spheres, nterm_A
select nterm_B, chain B and resi 1
color green, nterm_B
show spheres, nterm_B
select nterm_C, chain C and resi 1
color green, nterm_C
show spheres, nterm_C
select nterm_D, chain D and resi 1
color green, nterm_D
show spheres, nterm_D

# Show connections as distances
distance conn_A_D, nterm_A, nterm_D
distance conn_B_A, nterm_B, nterm_A
distance conn_D_C, nterm_D, nterm_C

# Color connected pairs in red
color red, chain A or chain D
color red, chain B or chain A
color red, chain D or chain C

# Summary: 3 connections from 4 chains
# Success rate: 75.0%

zoom all
