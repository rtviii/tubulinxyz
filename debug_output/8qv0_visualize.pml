# PyMOL script for 8QV0 N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/8QV0.cif)
load 8QV0.cif
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
select nterm_E, chain E and resi 1
color green, nterm_E
show spheres, nterm_E
select nterm_F, chain F and resi 1
color green, nterm_F
show spheres, nterm_F
select nterm_G, chain G and resi 1
color green, nterm_G
show spheres, nterm_G
select nterm_H, chain H and resi 1
color green, nterm_H
show spheres, nterm_H
select nterm_I, chain I and resi 1
color green, nterm_I
show spheres, nterm_I
select nterm_J, chain J and resi 1
color green, nterm_J
show spheres, nterm_J
select nterm_K, chain K and resi 1
color green, nterm_K
show spheres, nterm_K
select nterm_L, chain L and resi 1
color green, nterm_L
show spheres, nterm_L
select nterm_M, chain M and resi 1
color green, nterm_M
show spheres, nterm_M
select nterm_N, chain N and resi 1
color green, nterm_N
show spheres, nterm_N
select nterm_O, chain O and resi 1
color green, nterm_O
show spheres, nterm_O
select nterm_P, chain P and resi 1
color green, nterm_P
show spheres, nterm_P
select nterm_Q, chain Q and resi 1
color green, nterm_Q
show spheres, nterm_Q
select nterm_R, chain R and resi 1
color green, nterm_R
show spheres, nterm_R
select nterm_S, chain S and resi 1
color green, nterm_S
show spheres, nterm_S
select nterm_T, chain T and resi 1
color green, nterm_T
show spheres, nterm_T
select nterm_U, chain U and resi 1
color green, nterm_U
show spheres, nterm_U
select nterm_V, chain V and resi 1
color green, nterm_V
show spheres, nterm_V
select nterm_W, chain W and resi 1
color green, nterm_W
show spheres, nterm_W
select nterm_X, chain X and resi 1
color green, nterm_X
show spheres, nterm_X
select nterm_Y, chain Y and resi 1
color green, nterm_Y
show spheres, nterm_Y
select nterm_Z, chain Z and resi 1
color green, nterm_Z
show spheres, nterm_Z

# Show connections as distances
distance conn_B_A, nterm_B, nterm_A
distance conn_O_C, nterm_O, nterm_C
distance conn_P_D, nterm_P, nterm_D
distance conn_Q_E, nterm_Q, nterm_E
distance conn_R_F, nterm_R, nterm_F
distance conn_S_G, nterm_S, nterm_G
distance conn_T_H, nterm_T, nterm_H
distance conn_U_I, nterm_U, nterm_I
distance conn_V_J, nterm_V, nterm_J
distance conn_W_K, nterm_W, nterm_K
distance conn_X_L, nterm_X, nterm_L
distance conn_Y_M, nterm_Y, nterm_M
distance conn_Z_N, nterm_Z, nterm_N

# Color connected pairs in red
color red, chain B or chain A
color red, chain O or chain C
color red, chain P or chain D
color red, chain Q or chain E
color red, chain R or chain F
color red, chain S or chain G
color red, chain T or chain H
color red, chain U or chain I
color red, chain V or chain J
color red, chain W or chain K
color red, chain X or chain L
color red, chain Y or chain M
color red, chain Z or chain N

# Summary: 13 connections from 26 chains
# Success rate: 50.0%

zoom all
