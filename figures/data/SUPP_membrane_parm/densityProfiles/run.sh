#!/usr/bin/env bash

#FAITS: "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE:
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"
#traj456: "SCYM"
#(-1):"SCRN" "SCW"

aafile=("SCRN" "SCW")
profiles=("chains" "choline" "phosphate" "total" "water")
traj=(1 2 3)

for i in "${profiles[@]}"
do
  for aa in "${aafile[@]}"
  do
    for t in "${traj[@]}"
    do
      cp /media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-1/analyses/traj$t/data/densityProfiles/profile-${i}[3-5].dat .
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
    mv trajectory-dens.dat ${aa,,}-1-$i.dat
    if [ $aa = "NONE" ]; then
      mv ${aa,,}-$i.dat popc-$i.dat
    fi
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
