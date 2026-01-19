#!/bin/bash

DIR=../../charmm-gui/namd
PSFFILE=$DIR/step5_input.psf
DCDFILE=bilayer-centered.dcd

cp $DCDFILE ./centered.dcd
cp $PSFFILE ./bilayer.psf
vmd -dispdev text -e densityProfiles-dopc.vmd 
exit 0

