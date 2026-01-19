#! /bin/bash

# This script builds the PSF/PDB files of bilayer-solutes systems
# Up to 2 solutes are allowed in this script
# Other required scripts: packmol.sh, psfgen.inp and extraBondsFile.py
# Required files: original bilayer.psf/pdb files.



# Step 1: build side-chain analog systems with packmol
#   Distribute the solutes and SOD/CLA in a vacuum box, exclude bilayer center.
#   Infos to provide to script:
#	directory to find solute pdb files.
#	concentration of solutes.
#	"identity" of solutes (see pdb filenames from solute directory)
#	number of SOD/CLA from original bilayer system.
#	a/b/c from original bilayer system.
#	zexclude: membrane thickness to exclude from solute placement.
#	nwaters: number of waters from original bilayer system.
#	FOR CHARGED SOLUTES: adjust the number of CLA/SOD in script to neutralize systems.


#
# step 1: build side-chain analog systems with packmol
# Files created: solute pdb files, and sod.pdb and cla.pdb (required by step 2 below).
#
#mkdir ../../results/homoPOPC-aa/homoPOPC-sca


memb=MEMBRANE
sconc=0.2 #in M
sol1=SOLUTE
sol2="none"  # use "none" for no second solute
DIRESULT="../../results/${memb}-aa/${memb}-${sol1^^}"

mkdir -p $DIRESULT
rm -r $DIRESULT/*
cp ../../Lipid_Bilayer/${memb}/charmm-gui.tgz .
tar -zxvf charmm-gui.tgz
mv charmm-gui/charmm-gui-* ../
rm charmm-gui
mv charmm-gui-* charmm-gui

step7file="./charmm-gui/namd/step7_production.inp"



sed "s/sconc/$sconc/" packmol-${memb}.sh | sed "s/sol1/$sol1/" | sed "s/sol2/$sol2/"  > temp.sh
bash temp.sh

#
# step2 : combine solutes with bilayer system using VMD/psfgen
# Files created: bilayer.psf/pdb
#

##
## BUG:
## does not succeed at removing correctly the "#s2" from the file... makes script crash...
##
if [[ $sol2 != "none"  ]]; then
  sed "s/sol1/${sol1^^}/" psfgen.inp | sed "s/sol2/${sol2^^}/" | sed "s/\#s2//" > temp.inp
else
  sed "s/sol1/${sol1^^}/" psfgen.inp  > temp.inp
fi

#Put dcdfreq at 5000 instead of 50000
printf '%s\n' "./charmm-gui/namd/step7_production.inp"
# --------------------------------------------------------------------
# 1) Sauvegarde puis remplacement de dcdfreq
# --------------------------------------------------------------------
[[ -e "./charmm-gui/namd/step7_production.inp.bak" ]] || cp -- "$step7file" "${step7file}.bak"
sed -i -E 's/^(dcdfreq[[:space:]]+)[0-9]+;/\15000;/' "$step7file"


if [[ $sol1 == "scen" || $sol1 == "sche" || $sol1 == "schp" || $sol1 == "scp" ]]; then
  newline='parameters              toppar/top_all36_sidechains.str'
  # Ne rien faire si la ligne existe déjà
  if grep -q "^[[:space:]]*parameters[[:space:]]\+toppar/top_all36_sidechains\.str" "$step7file"; then
    echo "↩︎  « side‑chains » already in $step7file : no change done."
    return
  fi
  for i in 1 2 3 4 5 6; do
    step6file="./charmm-gui/namd/step6.${i}_equilibration.inp"
    last_param_lnum=$(grep -n -E '^[[:space:]]*parameters[[:space:]]' "$step6file" \
                    | tail -1 | cut -d: -f1)
    sed -i "${last_param_lnum}a $newline" "$step6file"
  done
  # Insérer juste avant la première directive « source »
  sed -i "/^[[:space:]]*source[[:space:]]/i $newline" "$step7file"
  echo "Line « side‑chains » added to $step7file"
  cp ../solutes/toppar/top_all36_sidechains.str ./charmm-gui/toppar/
  cp ../solutes/toppar-namd/top_all36_sidechains.str ./charmm-gui/namd/toppar/
  echo "Parameter file « side‑chains » added to $step7file"
else
  echo "sol1 = « $sol1 » : no parameter file added."
fi

vmd -dispdev text -e temp.inp

#
# step 3: build extrabonds file
# File created: extrabonds.txt (required by namd to avoid solute aggregation)

sed "s/sol1/${sol1^^}/" extraBondsFile.py | sed "s/sol2/${sol2^^}/" > temp.py
python temp.py bilayer-solutes.pdb > extrabonds.txt


#
# step4: remove temporary files
#

rm temp.inp temp.pdb temp.psf temp1.pdb temp1.psf temp.py temp.sh solute1.pdb solute2.pdb cla.pdb sod.pdb charmm-gui.tgz

#Need to fic GLPA for PMm
if [[ $memb == "PMm" ]]; then
  echo "Fixing GLPA for PMm membrane..."
  mv bilayer-solutes.pdb bugged.pdb
  bash ./fix_glpa.sh bugged.pdb charmm-gui/namd/step5_input.pdb > bilayer-solutes.pdb
  echo "Fixed GLPA for PMm membrane !"
fi

mv bilayer-solutes.pdb charmm-gui/namd/step5_input.pdb
mv bilayer-solutes.psf charmm-gui/namd/step5_input.psf
mv charmm-gui $DIRESULT/
cp ../../submit-* $DIRESULT/charmm-gui/namd/
echo "Done!"
