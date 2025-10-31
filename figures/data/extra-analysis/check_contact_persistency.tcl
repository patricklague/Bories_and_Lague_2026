# ===============================================================
# Script : contact_persistence.tcl
# Objectif : calculer la persistance des contacts entre analogues
#            "resname AAA" dans une trajectoire DCD
# Auteur : toi 😊 (avec un petit coup de main de ChatGPT)
# ===============================================================

# -----------------------
# Paramètres utilisateur
# -----------------------
set outfile [open "contact_persistence.dat" w]
set analog_sel "resname AAA"
set psffile "input.psf"
set trajfile "trajectory.dcd"
set cutoff 4.0

puts "=== DÉBUT DU SCRIPT ==="
puts "Chargement du fichier PSF : $psffile"
mol new $psffile type psf waitfor all

puts "Chargement de la trajectoire DCD : $trajfile"
mol addfile $trajfile type dcd waitfor all

# -----------------------
# Préparation des données
# -----------------------
set all [atomselect top "$analog_sel"]
set nframes [molinfo top get numframes]
set reslist [lsort -unique [$all get resid]]
set nres [llength $reslist]
puts "Nombre total de frames : $nframes"
puts "Nombre d'analogues trouvés : $nres"

# -----------------------
# Initialisation
# -----------------------
array set contact_duration {}
array set contact_count {}
array set durations_list {}
set total_contacts 0

puts "=== DÉBUT DE L'ANALYSE DES CONTACTS ==="

# -----------------------
# Boucle sur les frames
# -----------------------
for {set i 0} {$i < $nframes} {incr i} {
    molinfo top set frame $i
    puts "Analyse de la frame $i / $nframes ..."
    set contacts_this_frame {}
    set raw_contacts 0

    # --- Boucle sur toutes les paires d'analogues ---
    foreach r1 $reslist {
        set sel1 [atomselect top "$analog_sel and resid $r1"]
        foreach r2 $reslist {
            if {$r2 <= $r1} {continue}
            set sel2 [atomselect top "$analog_sel and resid $r2"]

            # Vérification du contact entre r1 et r2
            set pairs [measure contacts $cutoff $sel1 $sel2]
            set nAtomPairs [llength [lindex $pairs 0]]

            if {$nAtomPairs > 0} {
                lappend contacts_this_frame [list $r1 $r2]
                incr raw_contacts $nAtomPairs
            }

            $sel2 delete
        }
        $sel1 delete
    }

    # --- Affichage du nombre de contacts ---
    set nresPairs [llength $contacts_this_frame]
    puts "Frame $i : $nresPairs paires d'analogues en contact ; $raw_contacts paires d'atomes."
    incr total_contacts $nresPairs

    # --- Mise à jour des durées ---
    foreach pair $contacts_this_frame {
        lassign $pair r1 r2
        set key "$r1-$r2"
        if {[info exists contact_duration($key)]} {
            incr contact_duration($key)
        } else {
            set contact_duration($key) 1
        }
    }

    # --- Détection des contacts disparus ---
    foreach key [array names contact_duration] {
        set pair [split $key "-"]
        if {[lsearch -exact $contacts_this_frame $pair] == -1} {
            set d $contact_duration($key)
            if {$d > 1} {
                incr contact_count($key)
                lassign $pair r1 r2
                # 🔹 Enregistre la durée du contact pour chaque résidu
                lappend durations_list($r1) $d
                lappend durations_list($r2) $d
            }
            unset contact_duration($key)
        }
    }
}

# -----------------------
# Finalisation (contacts actifs à la fin)
# -----------------------
foreach key [array names contact_duration] {
    set d $contact_duration($key)
    if {$d > 1} {
        incr contact_count($key)
        lassign [split $key "-"] r1 r2
        lappend durations_list($r1) $d
        lappend durations_list($r2) $d
    }
}

puts "=== ANALYSE TERMINÉE ==="
puts "Écriture du fichier de sortie : contact_persistence.dat"

# -----------------------
# Écriture des résultats
# -----------------------
puts $outfile "#resid average_duration max_duration persistent_contacts"

foreach r $reslist {
    if {[info exists durations_list($r)]} {
        set durations $durations_list($r)
        set count [llength $durations]
        set sumdur 0
        foreach d $durations { incr sumdur $d }
        set avg [expr {$sumdur / double($count)}]
        set maxdur [tcl::mathfunc::max {*}$durations]
    } else {
        set avg 0
        set maxdur 0
        set count 0
    }
    puts $outfile "$r $avg $maxdur $count"
}

puts $outfile "#TOTAL_CONTACTS $total_contacts for $nframes frames"
close $outfile

puts "=== SCRIPT TERMINÉ AVEC SUCCÈS ==="
puts "Résultats écrits dans contact_persistence.dat"
puts "Total des contacts (toutes frames confondues) : $total_contacts"

exit


