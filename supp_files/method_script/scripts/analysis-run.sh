#!/bin/bash

# directory with the analysis scripts
#FAITS: "SCY" "SCF" "SCRN" "SCW" "SCI" "SCL" "SCV"
#A FAIRE: "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM"
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"

aafile=("SCY" "SCRN" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "GLYD" "SCV" "SCP")

nb_aa=26	#number of amino acids in the simulation (for contact analysis)

for aa in "${aafile[@]}"
do
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  traj=(1 2 3)
  if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
  elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
    traj=(4 2 3)
  fi
  for t in "${traj[@]}"
  do
    echo "Currently running $aa : trajectory $t" >> out.out
    cat analysis-self.sh | sed s=AMACID=$aa= | sed s=NUMTRAJ=$t= | sed s=NUMAA=$nb_aa= | sed s=RAWDATA=\"$DIR\"= > analysis-temp.sh
    bash analysis-temp.sh
    echo "$aa : traj $t  done" >> out.out
    rm analysis-temp.sh
  done
done
exit
