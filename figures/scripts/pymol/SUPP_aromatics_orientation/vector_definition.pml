cd Dropbox/SC_article/Bories_and_Lague_2025/figures/scripts/pymol/SUPP_aromatics_orientation/

load pdb/scf.pdb
load pdb/scy.pdb
#load pdb/scym.pdb
load pdb/scw.pdb
load pdb/sche.pdb
#load pdb/schd.pdb
#load pdb/schp.pdb

set ray_trace_mode, 3
show spheres
set sphere_scale, 0.3
set sphere_scale, 0.2, name H*
set stick_radius, 0.2

set sphere_color, gray50, name C*
set sphere_color, white, name H*
set sphere_color, blue, name N*
set sphere_color, red, name O*
set stick_color, white

set bg_rgb, white
select theta2_0, scf and (name CG or name CZ)
select theta2_1, scy and (name CG or name CZ)
select theta2_2, scw and (name CZ3 or name CE2)
select theta2_3, sche and (name CG or name CE1)
orient
zoom all, 0.5

set label_size, 36
set label_font_id, 5
#set label_color, back
set label_bg_color, back

cmd.label("theta2_0", "name")
cmd.label("theta2_1", "name")
cmd.label("theta2_2", "name")
cmd.label("theta2_3", "name")

#align scy and (name CG or name CB or name HB*), scf and (name CG or name CB or name HB*)
#align scw and (name CG or name CB or name HB*), scf and (name CG or name CB or name HB*)
#align sche and (name CG or name CB or name HB*), scf and (name CG or name CB or name HB*)

rotate z, 90
rotate y, 180
hide everything, scy or sche or scw
ray
png ../../../plot/SUPP_aromatics_orientation/illustration_scf_v2.png

hide everything, scf
show spheres, scy
show sticks, scy
cmd.label("all", "")
cmd.label("theta2_1", "name")
ray
png ../../../plot/SUPP_aromatics_orientation/illustration_scy_v2.png

hide everything, scy
show spheres, scw
show sticks, scw
cmd.label("all", "")
cmd.label("theta2_2", "name")
ray
png ../../../plot/SUPP_aromatics_orientation/illustration_scw_v2.png

hide everything, scw
show spheres, sche
show sticks, sche
cmd.label("all", "")
cmd.label("theta2_3", "name")
ray
png ../../../plot/SUPP_aromatics_orientation/illustration_sche_v2.png

#cmd.label("all", "")

rotate y, 90
ray
png ../../../plot/SUPP_aromatics_orientation/illustration_sche_v1.png

hide everything, sche
show spheres, scy
show sticks, scy
cmd.label("all", "")
cmd.label("theta2_1", "name")
ray
png ../../../plot/SUPP_aromatics_orientation/illustration_scy_v1.png

hide everything, scy
show spheres, scw
show sticks, scw
cmd.label("all", "")
cmd.label("theta2_2", "name")
ray
png ../../../plot/SUPP_aromatics_orientation/illustration_scw_v1.png

hide everything, scw
show spheres, scf
show sticks, scf
cmd.label("all", "")
cmd.label("theta2_0", "name")
ray
png ../../../plot/SUPP_aromatics_orientation/illustration_scf_v1.png




