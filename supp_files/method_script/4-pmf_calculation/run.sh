#!/usr/bin/env bash

#aafile=("SCA" "SCV" "SCL" "SCI") #hydrophobic side-chains
#aafile=("SCC" "SCM") #sulfured side-chains
#aafile=("SCS" "SCT" "SCN" "SCQ") #polar side-chains
#aafile=("SCF" "SCY" "SCW") #aromatic side-chains
#aafile=("PRO" "GLYD") #backbone-like side-chains
#aafile=("SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN") #titratable neutral side-chains
aafile=("SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR") #charged side-chains

for aa in "${aafile[@]}"
do
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  traj=(1 2 3)
  if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
  elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
    traj=(4 2 3)
  elif [[ "$aa" == "SCYM" ]]; then
    traj=(4 5 6)
  fi
  mkdir -p popc/distributions/${aa,,}
  mkdir -p popc/pmfs/
  mkdir -p pmf-data/
  for t in "${traj[@]}"
  do
    cp $DIR/analyses/traj$t/data/densityProfiles/profile-aa-600.dat ./trajectory$t.dat
    if [ $t = 4 ]; then
      mv trajectory$t.dat trajectory1.dat
    if [ $t = 5 ]; then
      mv trajectory$t.dat trajectory2.dat
    if [ $t = 6 ]; then
      mv trajectory$t.dat trajectory3.dat
    fi
  done
  python pmf-from-distribution-total.py
  mv trajectory*dat popc/distributions/${aa,,}/
  mv pmf_moyen.dat popc/pmfs/${aa,,}.dat
  mv pmf.dat pmf-data/pmf_${aa,,}.dat
done

