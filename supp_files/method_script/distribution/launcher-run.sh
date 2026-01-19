#!/bin/bash

# directory with the analysis scripts
#aafile=("GLYD" "SCP") #no orderparameters (traj1 GLYD déjà faite)
#aafile=("SCA" "SCV" "SCI" "SCL" "SCF" "SCY" "SCW") #no orderparameters
#aafile=("SCC" "SCS" "SCT") done # "SCM" "SCP" "SCHD" "SCHE" ) done #no orderparameters
#aafile=("SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN") #no orderparameters
#aafile=("SCDN" "SCEN" "SCQ" "SCN") #no orderparameters (traj1 GLYD déjà faite)
#aafile=("SCHP" "SCK" "SCR")	 #no orderparameters #group with -N at the end
#aafile=("SCCM" "SCD")	 #no orderparameters #group with -N, -05-N or -1-N at the end and at 50mM (6 aa) and 100mM (13 aa)
#aafile=("SCYM")
aafile=("GLYD" "SCP" "SCA" "SCL" "SCI" "SCHD" "SCHE" "SCHP" "SCF" "SCY" "SCYM" "SCS" "SCQ" "SCN" "SCM" "SCC" "SCCM" "SCK" "SCKN" "SCRN" "SCEN" "SCD" "SCDN" )
aafile=( "SCP" "SCHE" "SCHP" "SCEN")
aafile=("SCW" "SCRN")

#memb=("PMm" "G+PM" "G-IM" "homoPOPC")
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
