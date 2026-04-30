#!/usr/bin/env bash

aafile=( \
  "SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCN" "SCQ" \
  "SCF" "SCY" "SCW" "SCP" "GLYD" \
  "SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN" \
  "SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR" \
  "SCRN-1" "SCW-1" \
  "NONE" \
)

traj=(1 2 3)

# Canonical destination for the per-analog merged file
OUTDIR="../../../../../figures/data/SUPP_membrane_parm/area_per_lipid"
mkdir -p "$OUTDIR"

for aa in "${aafile[@]}"
do
  # Original (local) source of the per-trajectory data:
  DIR=../../../POPC-aa/POPC-$aa

  for t in "${traj[@]}"
  do
    cp $DIR/analyses/traj$t/data/cell.dat ./cell$t.dat
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

