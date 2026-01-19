# vmd -dispdev text -e pbc-recenter.tcl
package require pbctools

mol load psf temp.psf dcd temp.dcd

#pbc join res -ref "name CA"
#pbc wrap -centersel "segname PROA" -center com -compound res -all
#pbc wrap -centersel "segname MEMB and residue 1 to 200" -center com -compound res -all
#pbc wrap -centersel "segname MEMB" -center com -compound res -all

#pbc join res -ref "name CA"
#pbc wrap -centersel "segname PROA" -center com -compound res -all
#pbc wrap -centersel "segname MEMB and residue 1 to 200" -center com -compound res -all
pbc wrap -centersel "segname MEMB" -center com -compound res -all
#pbc wrap -centersel "segname PROA" -center com -compound res -all

# translate each frame to about Z=0
set total [atomselect top all]
#set ref [atomselect top "segname PROA or segname MEMB"]
set ref [atomselect top "segname MEMB"]
set nf [molinfo top get numframes]
puts [format "%i frames\n" $nf]
# Center the reference system around (0, 0, 0)
for {set f 0} {$f < $nf} {incr f} {
  $ref frame $f
  $total frame $f
  $total moveby [vecinvert [measure center $ref]]
}

animate write dcd centered.dcd

# save last frame to pdb
[atomselect top all frame last] writepdb lastFrame.pdb

exit
