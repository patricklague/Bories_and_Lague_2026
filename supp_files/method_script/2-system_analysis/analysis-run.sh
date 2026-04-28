#!/bin/bash

# directory with the analysis scripts
aafile=( \
  "SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCN" "SCQ" \
  "SCF" "SCY" "SCW" "SCP" "GLYD" \
  "SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN" \
  "SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR" \
  "SCRN-1" "SCW-1" \
  "NONE" \
)

for aa in "${aafile[@]}"
do
  #the $DIR variable refer here to the data produced by the simulations if launched using 1-system_generation
  DIR=../..
  #if using the concatened dcd from borealis, use the line below instead of the one above
  #DIR=../../../POPC-$aa
  traj=(1 2 3)
  for t in "${traj[@]}"
  do
    echo "Currently running $aa : trajectory $t" >> out.out
    cat analysis_per_system.sh | sed s=AMACID=$aa= | sed s=NUMTRAJ=$t= | sed s=RAWDATA=\"$DIR\"= > analysis-temp.sh
    #if using the concatened dcd from borealis, use the line below instead of the one above
    #cat analysis_borealis.sh | sed s=AMACID=$aa= | sed s=NUMTRAJ=$t= | sed s=RAWDATA=\"$DIR\"= > analysis-temp.sh
    bash analysis-temp.sh
    echo "$aa : traj $t  done" >> out.out
    rm analysis-temp.sh
  done
done
exit
