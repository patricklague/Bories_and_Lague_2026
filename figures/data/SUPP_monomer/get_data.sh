#!/bin/bash
# ===============================================================
# Script d'analyse par lots pour gérer la mémoire
# Traite 10 sections à la fois avec stride=1 (1000 frames/lot)
# Sections 400-1000 = 601 sections -> ~61 lots de 10 sections
# Total: 60100 frames par trajectoire
# ===============================================================

# Liste des analogues
# FAITS :  "SCY"
#aafile=("SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCQ" "SCN" "SCF" "SCW" "GLYD" "SCP" "SCCM" "SCYM" "SCE" "SCEN" "SCD" "SCDN" "SCK" "SCKN" "SCR" "SCRN" "SCHE" "SCHD" "SCHP")
aafile=("SCHP" "SCYM")

nb_aa=26

# Paramètres des lots
BATCH_SIZE=51       # Nombre de sections par lot
SECTION_START=400   # Première section
SECTION_END=1000    # Dernière section

for aa in "${aafile[@]}"
do
  # Déterminer le répertoire source
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  traj=(1 2 3)
  if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
  elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
    traj=(4 2 3)
  elif [[ "$aa" == "SCYM" ]]; then
    traj=(4 5 6)
  fi
  mkdir ${aa,,}
  
  for t in "${traj[@]}"
  do
    echo "=============================================="
    echo "Processing $aa : trajectory $t"
    echo "=============================================="
    if [ "$t" -eq 4 ]; then
      traj_num=1
    elif [ "$t" -eq 5 ]; then
      traj_num=2
    elif [ "$t" -eq 6 ]; then
      traj_num=3
    else
      traj_num="$t"
    fi
    mkdir ${aa,,}/traj${traj_num}

    # Numéro de trajectoire pour les noms de fichiers
    cp $DIR/analyses/traj${t}/data/com_classified_*.dat ./${aa,,}/traj${traj_num}/
  done
done

echo "=============================================="
echo "All analyses completed!"
echo "=============================================="I
