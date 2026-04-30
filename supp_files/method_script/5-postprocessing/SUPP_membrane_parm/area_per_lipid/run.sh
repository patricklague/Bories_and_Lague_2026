#!/usr/bin/env bash

#FAITS: "SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE: 
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"
#"NONE"

aafile=("SCYM")

traj=(4 5 6)

# Canonical destination for the per-analog merged file
OUTDIR="../../../../figures/data/SUPP_membrane_parm/area_per_lipid"
mkdir -p "$OUTDIR"

for aa in "${aafile[@]}"
do
  # Original (local) source of the per-trajectory data:
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  # Reviewer-facing alternative (uncomment if running from the public POPC-aa tree):
  #DIR="../../../../results/POPC-aa/POPC-$aa"

  for t in "${traj[@]}"
  do
    cp $DIR/analyses/traj$t/data/cell.dat ./cell$t.dat
    if [ $t = 4 ]; then
      mv cell$t.dat cell1.dat
    elif [ $t = 5 ]; then
      mv cell$t.dat cell2.dat
    elif [ $t = 6 ]; then
      mv cell$t.dat cell3.dat
    fi
  done
  cat get_apl.py | sed s=FILENAME="${aa,,}"= > get_temp.py
  python get_temp.py
  rm cell*.dat get_temp.py

  out="${aa,,}-apl.dat"
  if [ "$aa" = 'NONE' ]; then
    out="popc-apl.dat"
    mv none-apl.dat "$out"
  fi
  mv "$out" "$OUTDIR/$out"
done

