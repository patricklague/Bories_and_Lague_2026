# Script VMD Tcl pour analyser les analogues SC_ et les lipides POPC
# Utilisation : vmd -dispdev text -e analyze_monomers.tcl -args structure.psf trajectory.dcd

set psf analog-AAA.psf
set dcd trajectory.dcd
mol new $psf type psf waitfor all
mol addfile $dcd type dcd waitfor all

# Fichier de sortie
set outfile [open "monomers_percent.dat" w]
puts $outfile "#frame monomers_%" 
#popc_isolated_%"

# Sélections
set sel_SC_all [atomselect top "resname AAA and not hydrogen"]
#set sel_POPC_all [atomselect top "resname POPC"]

# --- Comptes & listes de resid (à partir de la frame 0)
animate goto 0
$sel_SC_all update
#$sel_POPC_all update

set SC_resids   [lsort -unique [$sel_SC_all get resid]]
#set POPC_resids [lsort -unique [$sel_POPC_all get resid]]

set n_SC   [llength $SC_resids]
#set n_POPC [llength $POPC_resids]

puts "Total SC_ residues: $n_SC"
#puts "Total POPC lipids: $n_POPC"

# Boucle sur toutes les frames
set nframes [molinfo top get numframes]
for {set i 0} {$i < $nframes} {incr i} {
    animate goto $i
    $sel_SC_all update
    #$sel_POPC_all update

    # ----- SC_ monomères (cutoff 8 Å)
    set monomer_count 0
    foreach resid_SC $SC_resids {
        # voisins SC_ à 8 Å du résidu courant
        set sel_neighbors [atomselect top "(resname AAA and not resid $resid_SC and not hydrogen) and within 10 of (resname AAA and resid $resid_SC and not hydrogen)"]
        if {[$sel_neighbors num] == 0} {
            incr monomer_count
        }
        $sel_neighbors delete
    }

    # ----- POPC isolés (cutoff 6 Å)
    #set popc_isolated_count 0
    #foreach resid_POPC $POPC_resids {
    #    set sel_neighbors [atomselect top "(resname POPC and not resid $resid_POPC) and within 6 of (resname POPC and resid $resid_POPC)"]
    #    if {[$sel_neighbors num] == 0} {
    #        incr popc_isolated_count
    #    }
    #    $sel_neighbors delete
    #}

    # --- Calcul des pourcentages ---
    set percent_SC [expr {100.0 * $monomer_count / $n_SC}]
    #set percent_POPC [expr {100.0 * $popc_isolated_count / $n_POPC}]

    # --- Écriture dans le fichier ---
    puts $outfile "$i $percent_SC" 
    #$percent_POPC"

    puts "Frame $i: AAA monomers=${percent_SC}%" 
    #| POPC isolated=${percent_POPC}%"
}

close $outfile
puts "Résultats enregistrés dans monomers_percent.dat"
exit


