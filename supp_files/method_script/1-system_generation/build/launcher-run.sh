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

memb=("homoPOPC")
for mm in "${memb[@]}"
do
  for aa in "${aafile[@]}"
  do
    echo "running in $mm : $aa" >> out.out
    cat run.sh | sed s=SOLUTE=\"${aa,,}\"= | sed s=MEMBRANE=\"$mm\"= > run-temp.sh
    bash run-temp.sh
    echo "$aa in $mm done" >> out.out
    rm run-temp.sh
  done
done
exit
