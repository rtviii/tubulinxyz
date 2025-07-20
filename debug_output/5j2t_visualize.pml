# PyMOL script for 5J2T N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/5J2T.cif)
load 5J2T.cif
hide everything
show cartoon
color grey, all

# Color N-terminus residues in green
select nterm_A, chain A and resi 1
color green, nterm_A
show spheres, nterm_A
select nterm_C, chain C and resi 1
color green, nterm_C
show spheres, nterm_C
select nterm_D, chain D and resi 1
color green, nterm_D
show spheres, nterm_D
select nterm_B, chain B and resi 1
color green, nterm_B
show spheres, nterm_B

# Show connections as distances
distance conn_C_B, nterm_C, nterm_B
distance conn_D_C, nterm_D, nterm_C
distance conn_B_A, nterm_B, nterm_A

# Color connected pairs in red
color red, chain C or chain B
color red, chain D or chain C
color red, chain B or chain A

# Summary: 3 connections from 4 chains
# Success rate: 75.0%

zoom all
