#!/usr/bin/env bash

#FAITS: "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE: 
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"
#traj456:"SCYM"

aafile=("SCRN")

traj=(1 2 3)


for aa in "${aafile[@]}"
do
  mkdir ${aa,,}
  for t in "${traj[@]}"
  do
    cp /media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa/charmm-gui/namd/out$t/section1000.dcd ./trajectory$t.dcd
    if [ $t = 4 ]; then
      mv trajectory$t.dcd trajectory1.dcd
    elif [ $t = 5 ]; then
      mv trajectory$t.dcd trajectory2.dcd
    elif [ $t = 6 ]; then
      mv trajectory$t.dcd trajectory3.dcd
    fi
  done
  for t in 1 2 3
  do
    cp -r trajectory$t.dcd traj.dcd
    cp ../input_systems/${aa,,}/step5_input.psf ./input.psf
    vmd -dispdev text -e get_frame.tcl
    mv last_frame.pdb end_traj$t.pdb
  done
  rm trajectory* input* traj.dcd
  mv end* ${aa,,}
done

