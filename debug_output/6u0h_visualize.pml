# PyMOL script for 6U0H N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/6U0H.cif)
load 6U0H.cif
hide everything
show cartoon
color grey, all

# Color tubulin chains
select alpha_chains, (
# Color N-terminus residues in green
select nterm_A, chain A and resi 1
color green, nterm_A
show spheres, nterm_A
select nterm_B, chain B and resi 1
color green, nterm_B
show spheres, nterm_B

# Show connections as distances
distance conn_B_A, nterm_B, nterm_A

# Color by connection strength
# Green = N-terminus residues
# Red = connected pairs
color red, chain B or chain A

# Summary: 1 connections from 2 chains
# Success rate: 50.0%

zoom all
