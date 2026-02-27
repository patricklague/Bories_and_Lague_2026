# Charger les fichiers
package require pbctools 
mol new input.psf
mol addfile traj.dcd waitfor all
#
# 1st, 2-step centering on membrane, first step around lipid #1, second step around whole membrane
#
pbc unwrap
pbc join connected -now
pbc wrap -centersel "segname MEMB and resid 1" -center com -compound residue -all 
pbc wrap -centersel "segname MEMB" -center com -compound residue -all 

#
# 2nd step: center lipids at (0,0,0)
#
set total [atomselect top all]
set ref [atomselect top "(segname MEMB)"]
set nf [molinfo top get numframes]
puts [format "%i frames\n" $nf]
# Center the reference system around (0, 0, 0)
for {set f 0} {$f < $nf} {incr f} {
  $ref frame $f
  $total frame $f
  $total moveby [vecinvert [measure center $ref weight mass]]
  set cell [pbc get -now]
  puts "Frame:$f $cell"
}
pbc unwrap
pbc wrap -centersel "segname MEMB" -center com -compound residue -all 
# Set the current frame to the last frame
animate goto [expr {$nf - 1}]

# Save the coordinates of the last frame to a PDB file (or other format)
set outputfile "last_frame.pdb"
set sel [atomselect top "all"] ;# Select all atoms in the top molecule
$sel writepdb $outputfile
#animate write dcd "traj_test.dcd"
$sel delete

puts "Dernière frame sauvegardée dans $outfile"

quit

