#!/usr/bin/env bash

#aafile=("SCA" "SCV" "SCL" "SCI") #hydrophobic side-chains
#aafile=("SCC" "SCM") #sulfured side-chains
#aafile=("SCS" "SCT" "SCN" "SCQ") #polar side-chains
#aafile=("SCF" "SCY" "SCW") #aromatic side-chains
#aafile=("SCP" "GLYD") #backbone-like side-chains
#aafile=("SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN") #titratable neutral side-chains
#aafile=("SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR") #charged side-chains

aafile=("NONE")
mode="monomer_4.5A"

for aa in "${aafile[@]}"
do
  DIR=../../../figures/data/distribution_data/$mode/${aa,,}
  traj=(1 2 3)
  
  #mkdir -p $output_dir/distributions_norm/${aa,,}
  mkdir -p $mode/pmfs/
  mkdir -p $mode/raw_data/
  cp $DIR/trajectory*.dat .

  python pmf-from-distribution.py
  #mv trajectory*dat popc/distributions/${aa,,}/
  rm trajectory*.dat
  mv pmf_moyen.dat $mode/pmfs/${aa,,}.dat
  mv pmf.dat $mode/raw_data/pmf_${aa,,}.dat
done

