#!/bin/bash
# ===============================================================
# Batch analysis script to manage memory
# Processes 10 sections at a time with stride=1 (1000 frames/batch)
# Sections 401-1000 = 600 sections -> 3 batches of 200 sections
# Total: 60000 frames per trajectory
# ===============================================================

# List of analogs
#aafile=("SCA" "SCV" "SCL" "SCI") #hydrophobic side-chains
#aafile=("SCC" "SCM") #sulfured side-chains
#aafile=("SCS" "SCT" "SCN" "SCQ") #polar side-chains
#aafile=("SCF" "SCY" "SCW") #aromatic side-chains
#aafile=("SCP" "GLYD") #backbone-like side-chains
#aafile=("SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN") #titratable neutral side-chains
#aafile=("SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR") #charged side-chains
aafile=("BBB")
traj=TTT
nb_aa=26

# Paramètres des lots
BATCH_SIZE=200       # Number of sections per batch
SECTION_START=401   # First section
SECTION_END=1000    # Last section

for aa in "${aafile[@]}"
do
  # Determine the source directory
  DIR=../../results/homoPOPC-aa/homoPOPC-$aa
  PSFFILE="$DIR/charmm-gui/namd/step5_input.psf"
  PDBFILE="$DIR/charmm-gui/namd/step5_input.pdb"
  # Copy the PSF file only once
  
  echo "=============================================="
  echo "Processing $aa : trajectory $traj"
  echo "=============================================="
  
    
  # Batch counting
  batch_num=0
    
  # Loop over sections in batches
  for (( start=$SECTION_START; start<=$SECTION_END; start+=$BATCH_SIZE ))
  do
    end=$((start + BATCH_SIZE - 1))
    if [ $end -gt $SECTION_END ]; then
      end=$SECTION_END
    fi
      
    echo "  Processing batch $batch_num: sections $start to $end"
      
    # Construct the list of DCD files for this batch
    files=""
    for (( i=$start; i<=$end; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$traj/section$i.dcd"
    done
    # Remove water
    cat ./indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp1.vmd
      
    vmd -dispdev text -e temp1.vmd  >& /dev/null
    PSFFILE2="./analog-$aa.psf"
    PDBFILE2="./analog-$aa.pdb"
    # Concatenate DCD files with stride=1 (all frames)
    catdcd -o trajectory.dcd -stride 1 -i findexfile.ind $files
      
    # Center the trajectory
    DCDFILE="trajectory.dcd"
    cat ./center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE=$DCDFILE= | sed s=AAA=$aa= > temp2.vmd
    vmd -dispdev text -e ./temp2.vmd >& /dev/null
    rm trajectory.dcd
    # Remove the membrane
    cat ./indexNoMemb.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=PDBFILE=$PDBFILE2= | sed s=AAA=$aa= > temp3.vmd
    vmd -dispdev text -e ./temp3.vmd >& /dev/null
    PSFFILE3="./analog-$aa.psf"
    catdcd -o trajectory.dcd -stride 1 -i findexfile.ind centered.dcd
    mv trajectory.dcd trajectory$traj.dcd
    rm trajectory.dcd centered.dcd temp*
  done
done

echo "=============================================="
echo "All analyses completed!"
echo "=============================================="I
