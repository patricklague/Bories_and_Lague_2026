#!/usr/bin/env bash

#FAITS: "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE:
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"
#traj456: "SCYM"
#(-1):"SCRN" "SCW"
#0.1M : "SCRN-1" "SCW-1"

aafile=("SCYM")
#profiles=("chains" "choline" "phosphate")
profiles=("total" "water")
traj=(1 2 3)

# Canonical destination for the per-analog merged file
OUTDIR="../../../../figures/data/SUPP_membrane_parm/densityProfiles"
mkdir -p "$OUTDIR"

for i in "${profiles[@]}"
do
  for aa in "${aafile[@]}"
  do
    # Original (local) source of the per-trajectory data:
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
    # Reviewer-facing alternative (uncomment if running from the public POPC-aa tree):
    #DIR="../../../../results/POPC-aa/POPC-$aa"

    for t in "${traj[@]}"
    do
      cp $DIR/analyses/traj$t/data/densityProfiles/profile-${i}[3-5].dat .
      for j in 3 4 5
      do
        mv profile-${i}$j.dat trajectory$t-$((j-2)).dat

        if [ $t = 4 ]; then
          mv trajectory$t-$((j-2)).dat trajectory1-$((j-2)).dat
        elif [ $t = 5 ]; then
          mv trajectory$t-$((j-2)).dat trajectory2-$((j-2)).dat
        elif [ $t = 6 ]; then
          mv trajectory$t-$((j-2)).dat trajectory3-$((j-2)).dat
        fi
      done
    done
    python get_density.py

    out="${aa,,}-$i.dat"
    if [ "$aa" = "NONE" ]; then
      out="popc-$i.dat"
    fi
    mv trajectory-dens.dat "$out"
    mv "$out" "$OUTDIR/$out"
    rm trajectory*.dat
  done
done

# When all densities have been calculated
#for i in "${profiles[@]}"
#do
#  cat calculate_mean_density.py | sed s=PROFILE=$i= > temp.py
#  python temp.py
#  mv all_dens.dat density-$i.dat
#done
