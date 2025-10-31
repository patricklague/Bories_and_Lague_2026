#!/bin/bash

# directory with the analysis scripts
#FAITS: "NONE"
#A FAIRE: "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM"
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" #dens_prof (en cours) à faire (SCE et SCD à faire)
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW" #dens_prof à faire
#aafile=("SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW") # a faire a partir de SCF
aafile=("SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ")
#aafile=("SCY")
nb_aa=26	#number of amino acids in the simulation (for contact analysis)

for aa in "${aafile[@]}"
do
  mkdir -p contacts/${aa,,}
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  traj=(1 2 3)
  if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
  elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
    traj=(4 2 3)
  fi
  for t in "${traj[@]}"
  do
    echo "Currently running $aa : trajectory $t"
    files=""
    for (( i=400; i<=1000; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    catdcd -o trajectory.dcd -stride 100 $files
    cp ../../../supp_files/input_systems/${aa,,}/step5_input.psf ./input.psf
    sed -r "s/AAA/$aa/g" get_purcentage.tcl > temp.tcl
    vmd -dispdev text -e temp.tcl
    sed -r "s/AAA/$aa/g" check_contact_persistency.tcl > temp.tcl
    vmd -dispdev text -e temp.tcl
    if [ $t = 4 ]; then
      mv monomers_percent.dat percent_traj1.dat
      mv contact_persistence.dat persistence_traj1.dat
    else
      mv monomers_percent.dat percent_traj$t.dat
      mv contact_persistence.dat persistence_traj$t.dat
    fi
    rm temp.tcl input.psf trajectory.dcd
  done
  bash summary_percent.sh
  bash summary_persistence.sh
  mv percent_summary.dat contacts/${aa,,}/percent_${aa,,}.dat
  mv summary_persistence.dat contacts/${aa,,}/persistence_${aa,,}.dat
  rm percent* contact_persistence*
done



