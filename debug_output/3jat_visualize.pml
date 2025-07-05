<<<<<<<< HEAD:debug_output/3jat_visualize.pml
# PyMOL script for 3JAT N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/3JAT.cif)
load 3JAT.cif
========
# PyMOL script for 6O2R N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/6O2R.cif)
load 6O2R.cif
>>>>>>>> origin/main:debug_output/6o2r_visualize.pml
hide everything
show cartoon
color grey, all

# Color N-terminus residues in green
select nterm_E, chain E and resi 1
color green, nterm_E
show spheres, nterm_E
select nterm_F, chain F and resi 1
color green, nterm_F
show spheres, nterm_F
select nterm_J, chain J and resi 1
color green, nterm_J
show spheres, nterm_J
select nterm_G, chain G and resi 1
color green, nterm_G
show spheres, nterm_G
select nterm_C, chain C and resi 1
color green, nterm_C
show spheres, nterm_C
select nterm_D, chain D and resi 1
color green, nterm_D
show spheres, nterm_D
select nterm_L, chain L and resi 1
color green, nterm_L
show spheres, nterm_L
select nterm_I, chain I and resi 1
color green, nterm_I
show spheres, nterm_I
select nterm_A, chain A and resi 1
color green, nterm_A
show spheres, nterm_A
select nterm_B, chain B and resi 1
color green, nterm_B
show spheres, nterm_B
select nterm_K, chain K and resi 1
color green, nterm_K
show spheres, nterm_K
select nterm_H, chain H and resi 1
color green, nterm_H
show spheres, nterm_H

# Show connections as distances
distance conn_E_F, nterm_E, nterm_F
distance conn_F_J, nterm_F, nterm_J
distance conn_G_E, nterm_G, nterm_E
distance conn_C_D, nterm_C, nterm_D
distance conn_D_L, nterm_D, nterm_L
distance conn_I_C, nterm_I, nterm_C
distance conn_A_B, nterm_A, nterm_B
distance conn_B_K, nterm_B, nterm_K
distance conn_H_A, nterm_H, nterm_A

# Color connected pairs in red
<<<<<<<< HEAD:debug_output/3jat_visualize.pml
========
color red, chain A or chain B
color red, chain B or chain K
color red, chain C or chain D
color red, chain D or chain L
>>>>>>>> origin/main:debug_output/6o2r_visualize.pml
color red, chain E or chain F
color red, chain F or chain J
color red, chain G or chain E
color red, chain C or chain D
color red, chain D or chain L
color red, chain I or chain C
color red, chain A or chain B
color red, chain B or chain K
color red, chain H or chain A

# Summary: 9 connections from 12 chains
# Success rate: 75.0%

zoom all
