#!/usr/bin/env bash

aafile=( \
  "SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCN" "SCQ" \
  "SCF" "SCY" "SCW" "SCP" "GLYD" \
  "SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN" \
  "SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR" \
  "SCRN-1" "SCW-1" \
  "NONE" \
)

# Canonical destination for the per-analog merged file
OUTDIR="../../../../../figures/data/SUPP_membrane_parm/thickness"
mkdir -p "$OUTDIR"

for aa in "${aafile[@]}"
do
  # Original (local) source of the per-trajectory data:
  DIR=../../../POPC-aa/POPC-$aa
  traj=(1 2 3)
  for t in "${traj[@]}"
  do
    cp $DIR/analyses/traj$t/data/thickness.dat ./thickness$t.dat
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

