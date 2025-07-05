# PyMOL script for 8VRK N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/8VRK.cif)
load 8VRK.cif
hide everything
show cartoon
color grey, all

# Color tubulin chains
select alpha_chains, (
# Color N-terminus residues in green
select nterm_1, chain 1 and resi 1
color green, nterm_1
show spheres, nterm_1
select nterm_2, chain 2 and resi 1
color green, nterm_2
show spheres, nterm_2
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
select nterm_o, chain o and resi 1
color green, nterm_o
show spheres, nterm_o
select nterm_p, chain p and resi 1
color green, nterm_p
show spheres, nterm_p
select nterm_q, chain q and resi 1
color green, nterm_q
show spheres, nterm_q
select nterm_r, chain r and resi 1
color green, nterm_r
show spheres, nterm_r
select nterm_s, chain s and resi 1
color green, nterm_s
show spheres, nterm_s
select nterm_t, chain t and resi 1
color green, nterm_t
show spheres, nterm_t
select nterm_u, chain u and resi 1
color green, nterm_u
show spheres, nterm_u
select nterm_v, chain v and resi 1
color green, nterm_v
show spheres, nterm_v
select nterm_w, chain w and resi 1
color green, nterm_w
show spheres, nterm_w
select nterm_x, chain x and resi 1
color green, nterm_x
show spheres, nterm_x
select nterm_y, chain y and resi 1
color green, nterm_y
show spheres, nterm_y
select nterm_z, chain z and resi 1
color green, nterm_z
show spheres, nterm_z

# Show connections as distances
distance conn_2_1, nterm_2, nterm_1
distance conn_o_O, nterm_o, nterm_O
distance conn_p_P, nterm_p, nterm_P
distance conn_q_Q, nterm_q, nterm_Q
distance conn_r_R, nterm_r, nterm_R
distance conn_s_S, nterm_s, nterm_S
distance conn_t_T, nterm_t, nterm_T
distance conn_u_U, nterm_u, nterm_U
distance conn_v_V, nterm_v, nterm_V
distance conn_w_W, nterm_w, nterm_W
distance conn_x_X, nterm_x, nterm_X
distance conn_y_Y, nterm_y, nterm_Y
distance conn_z_Z, nterm_z, nterm_Z

# Color by connection strength
# Green = N-terminus residues
# Red = connected pairs
color red, chain 2 or chain 1
color red, chain o or chain O
color red, chain p or chain P
color red, chain q or chain Q
color red, chain r or chain R
color red, chain s or chain S
color red, chain t or chain T
color red, chain u or chain U
color red, chain v or chain V
color red, chain w or chain W
color red, chain x or chain X
color red, chain y or chain Y
color red, chain z or chain Z

# Summary: 13 connections from 26 chains
# Success rate: 50.0%

zoom all
