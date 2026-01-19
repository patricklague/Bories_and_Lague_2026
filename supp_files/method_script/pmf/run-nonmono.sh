#!/usr/bin/env bash

#FAITS: "SCV" "GLYD" "SCA" "SCP" "SCW" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC"
#A FAIRE: "SCM"
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"

aafile=( "SCY" "SCF" "SCRN" "SCW" "SCI" "SCL") #SCV "SCY" "SCF" "SCRN" "SCW" "SCI" "SCL"

traj=(1 2 3)

aafile=("SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW")
for aa in "${aafile[@]}"
do
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  traj=(1 2 3)
  if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
  elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
    traj=(4 2 3)
  fi
  mkdir -p popc-nonmono/distributions/${aa,,}
  mkdir -p popc-nonmono/pmfs/
  mkdir -p pmf-data-nonmono/
  for t in "${traj[@]}"
  do
    cp $DIR/analyses/traj$t/data/densityProfiles/profile-aa-600-mono-2.dat ./trajectory$t.dat
    if [ $t = 4 ]; then
      mv trajectory$t.dat trajectory1.dat
    fi
  done
  python pmf-from-distribution6-nonmono.py
  mv trajectory*dat popc-nonmono/distributions/${aa,,}/
  mv pmf_moyen.dat popc-nonmono/pmfs/${aa,,}.dat
  mv pmf.dat pmf-data-nonmono/pmf_${aa,,}.dat
done

