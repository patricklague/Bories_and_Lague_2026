#!/usr/bin/env bash

#FAITS: "SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE: 
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"
#"NONE"

#aafile=("SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW")
aafile=("SCYM")



# Canonical destination for the per-analog merged file
OUTDIR="../../../../figures/data/SUPP_membrane_parm/thickness"
mkdir -p "$OUTDIR"

for aa in "${aafile[@]}"
do
  # Original (local) source of the per-trajectory data:
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  # Reviewer-facing alternative (uncomment if running from the public POPC-aa tree):
  #DIR="../../../../results/POPC-aa/POPC-$aa"
  traj=(1 2 3)
  if [[ "$aa" == "SCK" || "$aa" == "SCD" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
  elif [[ "$aa" == "SCYM" ]]; then
    traj=(4 5 6)
  fi
  for t in "${traj[@]}"
  do
    cp $DIR/analyses/traj$t/data/thickness.dat ./thickness$t.dat
    if [ $t = 4 ]; then
      mv thickness$t.dat thickness1.dat
    elif [ $t = 5 ]; then
      mv thickness$t.dat thickness2.dat
    elif [ $t = 6 ]; then
      mv thickness$t.dat thickness3.dat
    fi
  done
  cat get_thickness.py | sed s=FILENAME="${aa,,}"= > get_temp.py
  python get_temp.py
  rm thickness*.dat get_temp.py

  out="${aa,,}-thickness.dat"
  if [ "$aa" = 'NONE' ]; then
    out="popc-thickness.dat"
    mv none-thickness.dat "$out"
  fi
  mv "$out" "$OUTDIR/$out"
done

