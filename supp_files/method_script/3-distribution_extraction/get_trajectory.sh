#!/bin/bash
# ===============================================================
# Script d'analyse par lots pour gérer la mémoire
# Traite 10 sections à la fois avec stride=1 (1000 frames/lot)
# Sections 400-1000 = 601 sections -> ~61 lots de 10 sections
# Total: 60100 frames par trajectoire (au lieu de 601 avec stride=100)
# ===============================================================

# Liste des analogues
#aafile=("SCA" "SCV" "SCL" "SCI") #hydrophobic side-chains
#aafile=("SCC" "SCM") #sulfured side-chains
#aafile=("SCS" "SCT" "SCN" "SCQ") #polar side-chains
#aafile=("SCF" "SCY" "SCW") #aromatic side-chains
#aafile=("SCP" "GLYD") #backbone-like side-chains
#aafile=("SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN") #titratable neutral side-chains
#aafile=("SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR") #charged side-chains
aafile=("BBB")
nb_aa=26

# Paramètres des lots
BATCH_SIZE=600       # Nombre de sections par lot
SECTION_START=401   # Première section
SECTION_END=1000    # Dernière section

mkdir -p contacts/

for aa in "${aafile[@]}"
do
  # Déterminer le répertoire source
  DIR=../../results/POPC-aa/POPC-$aa
  traj=(1 2 3)
  PSFFILE="$DIR/charmm-gui/namd/step5_input.psf"
  PDBFILE="$DIR/charmm-gui/namd/step5_input.pdb"
  # Copy the PSF file only once
  
  for t in "${traj[@]}"
  do
    echo "=============================================="
    echo "Processing $aa : trajectory $t"
    echo "=============================================="
    
    # Trajectory number for file names
    DATA="$DIR/analyses/traj$t/data"
    
    # Batch counter
    batch_num=0
    
    # Loop over sections by batches
    for (( start=$SECTION_START; start<=$SECTION_END; start+=$BATCH_SIZE ))
    do
      end=$((start + BATCH_SIZE - 1))
      if [ $end -gt $SECTION_END ]; then
        end=$SECTION_END
      fi
      
      echo "  Processing batch $batch_num: sections $start to $end"
      
      # Build the list of DCD files for this batch
      files=""
      for (( i=$start; i<=$end; i++ ))
      do
        files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
      done
      #remove water
      cat ./indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp1.vmd
      
      vmd -dispdev text -e temp1.vmd  >& /dev/null
      PSFFILE2="./analog-$aa.psf"
      PDBFILE2="./analog-$aa.pdb"
      # Concatenate the DCD files with stride=1 (all frames)
      catdcd -o trajectory.dcd -stride 1 -i findexfile.ind $files
      
      #center the trajectory
      DCDFILE="trajectory.dcd"
      cat ./center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE=$DCDFILE= | sed s=AAA=$aa= > temp2.vmd
      vmd -dispdev text -e ./temp2.vmd >& /dev/null
      rm trajectory.dcd
      #remove the membrane
      cat ./indexNoMemb.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=PDBFILE=$PDBFILE2= | sed s=AAA=$aa= > temp3.vmd
      vmd -dispdev text -e ./temp3.vmd >& /dev/null
      PSFFILE3="./analog-$aa.psf"
      catdcd -o trajectory.dcd -stride 1 -i findexfile.ind centered.dcd
      mv trajectory.dcd trajectory$t.dcd
      rm trajectory.dcd centered.dcd temp*
    done
  done
done

echo "=============================================="
echo "All analyses completed!"
echo "=============================================="I
