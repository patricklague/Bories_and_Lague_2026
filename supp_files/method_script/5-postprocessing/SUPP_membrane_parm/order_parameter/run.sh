#!/usr/bin/env bash

#FAITS: "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE:
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"
#traj456: "SCYM"
#(-1):"SCRN" "SCW"

aafile=("SCD" "SCK")
chain=(2 3)
traj=(1 2 3)

# Canonical destination for the per-analog merged file
OUTDIR="../../../../../figures/data/SUPP_membrane_parm/order_parameter"
mkdir -p "$OUTDIR"

for i in "${chain[@]}"
do
  for aa in "${aafile[@]}"
  do
    # Original (local) source of the per-trajectory data:
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N

    for t in "${traj[@]}"
    do
      cp $DIR/analyses/traj$t/data/orderParameters/orderparameters-chain${i}[3-5].dat .
      for j in 3 4 5
      do
        mv orderparameters-chain${i}$j.dat scd${i}-t${t}-$((j-2)).dat
        if [ $t = 4 ]; then
          mv scd${i}-t${t}-$((j-2)).dat scd${i}-t1-$((j-2)).dat
        elif [ $t = 5 ]; then
          mv scd${i}-t${t}-$((j-2)).dat scd${i}-t2-$((j-2)).dat
        elif [ $t = 6 ]; then
          mv scd${i}-t${t}-$((j-2)).dat scd${i}-t3-$((j-2)).dat
        fi
      done
    done
    python get_scd.py

    if [ "$aa" = "NONE" ]; then
      out="popc-chain$i.dat"
    else
      out="${aa,,}-chain$i.dat"
    fi
    mv trajectory_scd.dat "$out"
    mv "$out" "$OUTDIR/$out"
    rm scd*-t*.dat
  done
done

