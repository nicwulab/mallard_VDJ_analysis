#!/usr/bin/env pymol
load result/TMalign/fold_rag1_rag2_human_model_1.pdb_all_atm, format=pdb
hide all
show cartoon
color blue, chain A
color red, chain B
set ray_shadow, 0
set stick_radius, 0.3
set sphere_scale, 0.25
show stick, not polymer
show sphere, not polymer
bg_color white
set transparency=0.2
zoom polymer

