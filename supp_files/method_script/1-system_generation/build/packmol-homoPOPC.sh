#!/bin/bash


#######################################
# number of solute molecules
#######################################
#
# Sc: solute concentration (mol/L)
#
# Other parameters are from file step5_assembly.str
# A, B, and C (box dimensions)
# and number of water molecules
# zexclude= +/-zdim to exclude solute (bilayer slab)

pdbdir="../solutes/"
Sc=sconc
solute1="sol1"
solute2="sol2"
sod="18"
cla="18"

# edit only if changing lipid bilayer
a="46.7504011" 
b="46.7504011"
c="145"
zexclude="22"
nwaters="7080"

#
# nsolute = (18.01528g/mol) x (1cc/0,99567g) x (1L/1000 cc) x (Sc mol/L) x (nwaters) 
# nsolute = 0,018093625 x (Sc mol/L) x (nwaters)
ns=$(echo "0.018093625 * $Sc * $nwaters" | bc)
nsolute=$(echo "($ns+0.5)/1" | bc)

# comment/uncomment to get neutralized systems
if [[ $solute1 == "scym" || $solute1 == "sccm" || $solute1 == "sce" || $solute1 == "scd" ]]; then
  nsod=$(echo "($nsolute+$sod)" | bc) #for negativ residue
else
  nsod=$(echo "($sod)" | bc)
fi

if [[ $solute1 == "schp" || $solute1 == "sck" || $solute1 == "scr" ]]; then
  ncla=$(echo "($nsolute+$cla)" | bc) #for positiv residue
else
  ncla=$(echo "($cla)" | bc)
fi

#nsod=$(echo "($sod)" | bc)
#nsod=$(echo "($nsolute+$sod)" | bc) #uncomment if negativ residue

#ncla=$(echo "($nsolute+$cla)" | bc) #uncomment if positiv residue
#ncla=$(echo "($cla)" | bc)


#### do not edit below #####
dimx=$(echo "$a/2.0" | bc)
dimy=$(echo "$b/2.0" | bc)
dimz=$(echo "$c/2.0" | bc)

#cat packmol-template.inp > temp.inp
rm -f temp.inp

echo "" >> temp.inp
echo "tolerance 2.0" >> temp.inp
echo "seed -1" >> temp.inp
echo "filetype pdb" >> temp.inp
echo "output packmol.pdb" >> temp.inp


if [[ $solute1 != "none"  ]]; then
  echo "" >> temp.inp
  echo "structure $pdbdir/$solute1.pdb" >> temp.inp
  echo "  number $nsolute " >> temp.inp
  echo "  inside box -$dimx -$dimy -$dimz $dimx $dimy $dimz" >> temp.inp
  echo "  outside box -$dimx -$dimy -$zexclude $dimx $dimy $zexclude" >> temp.inp
  echo "end structure"  >> temp.inp
  echo "" >> temp.inp
fi

#if [[ -v solute2  ]]; then
if [[ $solute2 != "none"  ]]; then
  echo "" >> temp.inp
  echo "structure $pdbdir/$solute2.pdb" >> temp.inp
  echo "  number $nsolute " >> temp.inp
  echo "  inside box -$dimx -$dimy -$dimz $dimx $dimy $dimz" >> temp.inp
  echo "  outside box -$dimx -$dimy -$zexclude $dimx $dimy $zexclude" >> temp.inp
  echo "end structure"  >> temp.inp
  echo "" >> temp.inp
fi

echo "" >> temp.inp
echo "structure $pdbdir/sod.pdb" >> temp.inp
echo "  number $nsod " >> temp.inp
echo "  inside box -$dimx -$dimy -$dimz $dimx $dimy $dimz" >> temp.inp
echo "  outside box -$dimx -$dimy -$zexclude $dimx $dimy $zexclude" >> temp.inp
echo "end structure"  >> temp.inp
echo "" >> temp.inp

echo "" >> temp.inp
echo "structure $pdbdir/cla.pdb" >> temp.inp
echo "  number $ncla " >> temp.inp
echo "  inside box -$dimx -$dimy -$dimz $dimx $dimy $dimz" >> temp.inp
echo "  outside box -$dimx -$dimy -$zexclude $dimx $dimy $zexclude" >> temp.inp
echo "end structure"  >> temp.inp
echo "" >> temp.inp

packmol < temp.inp > packmol.log

#grep SCF packmol.pdb > scf.pdb
grep ${solute1^^} packmol.pdb > solute1.pdb

#if [[ -v solute2  ]]; then
if [[ $solute2 != "none"  ]]; then
  grep ${solute2^^} packmol.pdb > solute2.pdb
fi

grep SOD packmol.pdb > sod.pdb
grep CLA packmol.pdb > cla.pdb





