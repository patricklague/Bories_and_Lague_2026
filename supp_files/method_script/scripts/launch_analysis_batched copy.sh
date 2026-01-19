#!/bin/bash
# ===============================================================
# Script d'analyse par lots pour gérer la mémoire
# Traite 10 sections à la fois avec stride=1 (1000 frames/lot)
# Sections 400-1000 = 601 sections -> ~61 lots de 10 sections
# Total: 60100 frames par trajectoire (au lieu de 601 avec stride=100)
# ===============================================================

# Liste des analogues
# FAITS :  "SCY"
aafile=("SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCQ" "SCN" "SCF" "SCW" "GLYD" "SCP" "SCCM" "SCYM" "SCE" "SCEN" "SCD" "SCDN" "SCK" "SCKN" "SCR" "SCRN" "SCHE" "SCHD" "SCHP")
#aafile=("SCI")  # Pour tester avec un seul analogue
nb_aa=26

# Paramètres des lots
BATCH_SIZE=50       # Nombre de sections par lot
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
  fi
  PSFFILE="$DIR/charmm-gui/namd/step5_input.psf"
  PDBFILE="$DIR/charmm-gui/namd/step5_input.pdb"
  # Copier le fichier PSF une seule fois
  cp ../../../supp_files/input_systems/${aa,,}/step5_input.psf ./input.psf
  
  for t in "${traj[@]}"
  do
    echo "=============================================="
    echo "Processing $aa : trajectory $t"
    echo "=============================================="
    
    # Numéro de trajectoire pour les noms de fichiers
    DATA="$DIR/analyses/traj$t/data"
    
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
      #retirer l'eau
      cat ./indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp1.vmd
      
      vmd -dispdev text -e temp1.vmd  >& /dev/null
      PSFFILE2="./analog-$aa.psf"
      PDBFILE2="./analog-$aa.pdb"
      # Concaténer les DCD avec stride=1 (toutes les frames)
      catdcd -o trajectory.dcd -stride 1 -i findexfile.ind $files
      
      #centrer la trajectoire
      DCDFILE="trajectory.dcd"
      cat ./center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE=$DCDFILE= | sed s=AAA=$aa= > temp2.vmd
      vmd -dispdev text -e ./temp2.vmd >& /dev/null
      rm trajectory.dcd
      #supprimer la membrane
      cat ./indexNoMemb.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=PDBFILE=$PDBFILE2= | sed s=AAA=$aa= > temp3.vmd
      vmd -dispdev text -e ./temp3.vmd >& /dev/null
      PSFFILE3="./analog-$aa.psf"
      catdcd -o trajectory.dcd -stride 1 -i findexfile.ind centered.dcd
      
      # Préparer le script TCL
      cat ./densityProfiles-aa.vmd | sed s=AAA=$aa=  | sed s=PSFFILE=$PSFFILE3= | sed s=DCDFILE="./trajectory.dcd"= > solutes.vmd
      
      # Exécuter l'analyse VMD
      vmd -dispdev text -e solutes.vmd
      # Sauvegarder les résultats du lot avec le numéro de frame global
      # Le numéro de frame global = batch_num * BATCH_SIZE * 100 + frame_local
      # (100 frames par section DCD)
      python ./average_profiles.py dens_mono_frame_*.dat > dens_mono_traj_${t}_batch_${batch_num}.dat
      python ./average_profiles.py dens_nonmono_frame_*.dat > dens_nonmono_traj_${t}_batch_${batch_num}.dat
      rm dens_mono_frame_*.dat dens_nonmono_frame_*.dat
      
      # Nettoyer la trajectoire temporaire
      rm -f trajectory.dcd solutes.vmd centered.dcd analog*
      
      batch_num=$((batch_num + 1))
    done
    python ./average_profiles.py dens_mono_traj_${t}_batch*.dat > dens_mono_avg.dat
    python ./average_profiles.py dens_nonmono_traj_${t}_batch*.dat > dens_nonmono_avg.dat
    ./combine_profiles.sh dens_mono_avg.dat dens_nonmono_avg.dat > dens_all_avg.dat
    mv -f dens_all_avg.dat  $DATA/densityProfiles/profile-aa-600-mono-2.dat
    mkdir contacts/${aa,,}
    mv dens_*mono_traj_${t}_batch*.dat contacts/${aa,,}/
  done
done

echo "=============================================="
echo "All analyses completed!"
echo "=============================================="I
