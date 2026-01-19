#!/bin/bash

###################################
#
# Script d'analyse membranaire
#
###################################

dir=../../charmm-gui/namd
PSFFILE=$dir/step5_input.psf
PDBFILE=$dir/bilayer.pdb
SCRIPTDIR=./
WORDOM=/usr/local/bin/wordom
CATDCD=/usr/local/bin/catdcd

#######################################
# what to do...
#######################################

# build the trajectory for the analyses
buildDcd="1"
center="1"      # center the peptide/bilayer in middle of system (required for almost all analyses)
cleanDcd="0"    # keep only the latest DCD file created, remove the previous ones.
first=1         # first section to include in analyses
last=190        # last section
stride=10       # make sure enough data for stats
                # usually stride=100 for whole trajectory
                # and stride=10-50 for production...

thickness="1"   # thickness of bilayer timeseries



########################
# do not edit below
########################

mkdir -p results

#######################################
# Build dcd trajectory without waters
#######################################

if [ "$buildDcd" = "1" ]; then

  echo ""
  echo "Building DCD trajectory from sections..."
  echo "Final files of this step: bilayer.pdb/psf/dcd"
  files=""
  for i in `seq 1 1`;
    do
      COUNTER=$first
      while [  $COUNTER -le $last ]; do
        if test -f "$dir/out/section$COUNTER.dcd"; then
          files="$files $dir/out/section$COUNTER.dcd"
        else
          echo "File $dir/out/section$COUNTER.dcd does not exists... exit!"
          exit
        fi
        echo $i $COUNTER
        let COUNTER=COUNTER+1
      done
    done

  $CATDCD -o bilayer.dcd -stride $stride  $files  >& /dev/null

  echo "     done!"
  echo ""
fi


#####################
# center the peptide
#####################

if [ "$center" = "1" ]; then

  echo ""
  echo "Centering of the bilayer ... bilayer-centered.dcd"

  if test -f bilayer.dcd ; then
	cp $PSFFILE ./bilayer.psf
    ln -s bilayer.psf temp.psf  # cheap friday PM trick...
    ln -s bilayer.dcd temp.dcd
    vmd -dispdev text -e $SCRIPTDIR/pbc-recenter.tcl >& /dev/null
    cp centered.dcd bilayer-centered.dcd
    rm temp.psf temp.dcd

    if [ "$cleanDcd" = "1" ]; then
      rm bilayer.dcd
    fi

    echo "Writing PDB of the last frame of the trajectory... results/bilayer-final.pdb"
    mv lastFrame.pdb results/bilayer-final.pdb

  else
    echo "There is no trajectory file to recenter!"
    echo ""
  fi

  echo "     done!"
  echo ""

fi

###############################
# Bilayer thickness timeseries
###############################

if [ "$thickness" = "1" ]; then

  echo ""
  echo "Bilayer thickness time series... bilayer should be centered before (to avoid splitting across PBC)!"
  echo "Calculated using VMD Membrane Plugin, with Step=10 and -sel {name P}"
  echo ""
  echo "Thickness values are stored in bilayer-thickness.dat"
  echo ""

  if test -f bilayer-centered.dcd ; then
    vmd -dispdev text -e $SCRIPTDIR/thickness.vmd
    # correct the header of results file...
    line=`head -n1 bilayer-thickness.dat`
    echo "#$line" > results/bilayer-thickness.dat
    tail -n+2 bilayer-thickness.dat >> results/bilayer-thickness.dat
    rm bilayer-thickness.dat
  else
    echo "There is no trajectory file to analyse!"
    echo ""
  fi

  echo "     done!"
  echo ""

fi
