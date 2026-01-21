#!/bin/bash
# ===============================================================
# Script d'analyse par lots pour gérer la mémoire
# Traite 10 sections à la fois avec stride=1 (1000 frames/lot)
# Sections 400-1000 = 601 sections -> ~61 lots de 10 sections
# Total: 60100 frames par trajectoire (au lieu de 601 avec stride=100)
# ===============================================================

# Liste des analogues
aafile=("SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW")
#aafile=("SCI")  # Pour tester avec un seul analogue
nb_aa=26

# Paramètres des lots
BATCH_SIZE=10       # Nombre de sections par lot
SECTION_START=400   # Première section
SECTION_END=1000    # Dernière section

for aa in "${aafile[@]}"
do
  mkdir -p contacts/${aa,,}
  mkdir -p contacts/${aa,,}/batches  # Dossier pour fichiers intermédiaires
  
  # Déterminer le répertoire source
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  traj=(1 2 3)
  if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
  elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
    traj=(4 2 3)
  fi
  
  # Copier le fichier PSF une seule fois
  cp ../../../supp_files/input_systems/${aa,,}/step5_input.psf ./input.psf
  
  for t in "${traj[@]}"
  do
    echo "=============================================="
    echo "Processing $aa : trajectory $t"
    echo "=============================================="
    
    # Numéro de trajectoire pour les noms de fichiers
    if [ $t = 4 ]; then
      traj_num=1
    else
      traj_num=$t
    fi
    
    # Compteur de lot
    batch_num=0
    
    # Boucle sur les sections par lots
    for (( start=$SECTION_START; start<=$SECTION_END; start+=$BATCH_SIZE ))
    do
      end=$((start + BATCH_SIZE - 1))
      if [ $end -gt $SECTION_END ]; then
        end=$SECTION_END
      fi
      
      echo "  Processing batch $batch_num: sections $start to $end"
      
      # Construire la liste des fichiers DCD pour ce lot
      files=""
      for (( i=$start; i<=$end; i++ ))
      do
        files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
      done
      
      # Concaténer les DCD avec stride=1 (toutes les frames)
      catdcd -o trajectory.dcd -stride 1 $files
      
      # Préparer le script TCL
      sed -r "s/AAA/$aa/g" get_purcentage.tcl > temp.tcl
      
      # Exécuter l'analyse VMD
      vmd -dispdev text -e temp.tcl
      
      # Sauvegarder les résultats du lot avec le numéro de frame global
      # Le numéro de frame global = batch_num * BATCH_SIZE * 100 + frame_local
      # (100 frames par section DCD)
      mv monomers_percent.dat contacts/${aa,,}/batches/percent_traj${traj_num}_batch${batch_num}.dat
      
      # Analyse de persistance des contacts
      sed -r "s/AAA/$aa/g" check_contact_persistency.tcl > temp.tcl
      vmd -dispdev text -e temp.tcl
      mv contact_persistence.dat contacts/${aa,,}/batches/persistence_traj${traj_num}_batch${batch_num}.dat
      
      # Nettoyer la trajectoire temporaire
      rm -f trajectory.dcd temp.tcl
      
      batch_num=$((batch_num + 1))
    done
    
    echo "  Merging batch files for trajectory $traj_num..."
    
    # Fusionner tous les lots pour cette trajectoire
    bash merge_batches.sh ${aa,,} $traj_num $batch_num
    
  done
  
  rm -f input.psf
  
  # Synthèse des 3 trajectoires
  echo "Computing final summary for $aa..."
  cd contacts/${aa,,}
  bash ../../summary_percent_full.sh
  bash ../../summary_persistence_full.sh
  cd ../..
  
  mv contacts/${aa,,}/percent_summary.dat contacts/${aa,,}/percent_${aa,,}.dat
  mv contacts/${aa,,}/persistence_summary.dat contacts/${aa,,}/persistence_${aa,,}.dat
  
done

echo "=============================================="
echo "All analyses completed!"
echo "=============================================="
