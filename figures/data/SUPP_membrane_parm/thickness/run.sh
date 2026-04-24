#!/usr/bin/env bash

#FAITS: "SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE: 
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"
#"NONE"

#aafile=("SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW")
aafile=("SCRN")

#traj=(4 5 6)


for aa in "${aafile[@]}"
do
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  traj=(1 2 3)
#  if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
#    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
#  elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
#    traj=(4 2 3)
#  elif [[ "$aa" == "SCYM" ]]; then
#    traj=(4 5 6)
#  fi
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
  rm thickness*.dat
  if [ $aa = 'NONE' ]; then
    mv none-thickness.dat popc-thickness.dat
  fi
done

