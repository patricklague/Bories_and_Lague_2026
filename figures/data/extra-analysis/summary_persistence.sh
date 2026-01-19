#!/bin/bash
# ===============================================================
# Script : summary_contact_persistence.sh
# Combine les fichiers contact_persistence.dat de plusieurs trajs
# et conserve les lignes #TOTAL_CONTACTS
# Sortie : summary_contact_persistence.dat
# ===============================================================

# Liste des fichiers à fusionner
files=("persistence_traj1.dat" "persistence_traj2.dat" "persistence_traj3.dat")

# Fichier de sortie
outfile="summary_persistence.dat"

# Vérification de l'existence des fichiers
for f in "${files[@]}"; do
    if [ ! -f "$f" ]; then
        echo "❌ Fichier $f introuvable, arrêt."
        exit 1
    fi
done

# En-tête du fichier résumé
echo "#resid avg_traj1 avg_traj2 avg_traj3 max_traj1 max_traj2 max_traj3 pers_traj1 pers_traj2 pers_traj3" > "$outfile"

# Extraction des lignes de données (sans les #TOTAL) et fusion
paste "${files[@]}" | awk '
    # On saute les lignes #TOTAL_CONTACTS
    $0 !~ /^#/ {
        resid=$1
        avg1=$2; max1=$3; pers1=$4
        avg2=$6; max2=$7; pers2=$8
        avg3=$10; max3=$11; pers3=$12
        printf "%d %.3f %.3f %.3f %d %d %d %d %d %d\n", resid, avg1, avg2, avg3, max1, max2, max3, pers1, pers2, pers3
    }
' >> "$outfile"

echo "#TOTAL traj1 traj2 traj3 nframes" >> "$outfile"
# ==============================================
# Traitement des lignes #TOTAL_CONTACTS
# ==============================================
total1=$(grep "#TOTAL_CONTACTS" "${files[0]}" | awk '{print $2}')
total2=$(grep "#TOTAL_CONTACTS" "${files[1]}" | awk '{print $2}')
total3=$(grep "#TOTAL_CONTACTS" "${files[2]}" | awk '{print $2}')
nframes=$(grep "#TOTAL_CONTACTS" "${files[0]}" | awk '{print $4}')

# Si un total est vide, on le remplace par 0
total1=${total1:-0}
total2=${total2:-0}
total3=${total3:-0}

# Écriture du résumé global
echo "#TOTAL $total1 $total2 $total3 $nframes" >> "$outfile"

echo "✅ Résumé global enregistré dans $outfile"

