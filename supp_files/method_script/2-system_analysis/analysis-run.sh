#!/bin/bash

# directory with the analysis scripts
#aafile=("SCA" "SCV" "SCL" "SCI") #hydrophobic side-chains
#aafile=("SCC" "SCM") #sulfured side-chains
#aafile=("SCS" "SCT" "SCN" "SCQ") #polar side-chains
#aafile=("SCF" "SCY" "SCW") #aromatic side-chains
#aafile=("PRO" "GLYD") #backbone-like side-chains
#aafile=("SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN") #titratable neutral side-chains
#aafile=("SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR") #charged side-chains
aafile=("NONE") #use "NONE" for the membrane only system



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
  for t in "${traj[@]}"
  do
    echo "Currently running $aa : trajectory $t" >> out.out
    cat analysis-self.sh | sed s=AMACID=$aa= | sed s=NUMTRAJ=$t= | sed s=RAWDATA=\"$DIR\"= > analysis-temp.sh
    bash analysis-temp.sh
    echo "$aa : traj $t  done" >> out.out
    rm analysis-temp.sh
  done
done
exit
