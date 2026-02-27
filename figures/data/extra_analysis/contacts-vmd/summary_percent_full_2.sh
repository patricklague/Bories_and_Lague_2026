#!/bin/bash
# ===============================================================
# Méthode rigoureuse : 1 trajectoire = 1 réplicat indépendant
# Calcule le % global par trajectoire, puis moyenne ± SE (n=3)
# 
# Cette méthode évite le problème d'autocorrélation des frames
# Utilise les fichiers bruts : contacts/{sc}/batches/percent_traj*_batch*.dat
# Format des fichiers : #frame monomers_% popc_isolated_%
# ===============================================================

echo "=== Méthode rigoureuse : 1 trajectoire = 1 réplicat ==="
echo ""

# Fichier de sortie
OUTPUT="percent_summary_rigorous.dat"

# En-tête du fichier de sortie
echo "# Méthode rigoureuse : moyenne ± erreur standard sur 3 réplicats indépendants" > $OUTPUT
echo "# Chaque trajectoire est traitée comme UNE mesure indépendante" >> $OUTPUT
echo "# analog monomer_mean monomer_se n_traj n_frames_total" >> $OUTPUT

# Liste des analogues (à adapter selon vos données)
ANALOGS="glyd scp sca scv scl sci scc scm scf scy scw scs sct scn scq scym sccm sce scen scd scdn scr scrn sck sckn sche schd schp"

for sc in $ANALOGS; do
    SC_UPPER=$(echo $sc | tr '[:lower:]' '[:upper:]')
    
    # Vérifier si le dossier batches existe
    if [ ! -d "contacts/$sc/batches" ]; then
        echo "  $SC_UPPER: dossier contacts/$sc/batches non trouvé, ignoré"
        continue
    fi
    
    echo "Traitement de $SC_UPPER..."
    
    # Initialiser les tableaux pour stocker les % par trajectoire
    mono_vals=""
    n_traj=0
    total_frames_all=0
    
    for traj in 1 2 3; do
        # Collecter tous les fichiers batch pour cette trajectoire
        FILES=$(ls contacts/$sc/batches/percent_traj${traj}_batch*.dat 2>/dev/null)
        #ls $FILES 
        if [ -z "$FILES" ]; then
            echo "  traj$traj: aucun fichier trouvé"
            continue
        fi
        
        # Calculer la moyenne des pourcentages de monomères sur tous les batches/frames
        # En accumulant toutes les valeurs puis en faisant la moyenne
        sum_mono=0
        n_frames=0
        
        for FILE in $FILES; do
            # Lire les pourcentages de monomères (colonne 2) en ignorant l'en-tête
            while read -r line; do
                # Ignorer les lignes de commentaire
                [[ "$line" =~ ^#.* ]] && continue
                [[ -z "$line" ]] && continue
                
                # Extraire le pourcentage de monomères (colonne 2)
                pct=$(echo "$line" | awk '{print $2}')
                
                if [ -n "$pct" ]; then
                    sum_mono=$(echo "$sum_mono + $pct" | bc -l)
                    n_frames=$((n_frames + 1))
                fi
            done < "$FILE"
        done
        
        # Si on a des données pour cette trajectoire
        if [ $n_frames -gt 0 ]; then
            # Calculer le pourcentage moyen pour cette trajectoire
            pct_mono=$(echo "scale=6; $sum_mono / $n_frames" | bc)
            
            mono_vals="$mono_vals $pct_mono"
            n_traj=$((n_traj + 1))
            total_frames_all=$((total_frames_all + n_frames))
            
            echo "  traj$traj: mono_mean=${pct_mono}% (n_frames=$n_frames)"
        fi
    done
    
    # Calculer moyenne et erreur standard si on a au moins 2 trajectoires
    if [ $n_traj -ge 2 ]; then
        # Utiliser awk pour les calculs statistiques
        result=$(echo "$mono_vals | $n_traj | $total_frames_all" | awk -F'|' '{
            # Parser les valeurs
            n = $2 + 0
            total_frames = $3 + 0
            
            # Monomères
            split($1, m, " ")
            sum_m = 0; sum_m2 = 0
            for (i in m) { if (m[i] != "") { sum_m += m[i]; sum_m2 += m[i]^2 } }
            mean_m = sum_m / n
            var_m = (sum_m2 - n * mean_m^2) / (n - 1)
            se_m = sqrt(var_m / n)
            
            printf "%.2f %.2f %d %d", mean_m, se_m, n, total_frames
        }')
        
        echo "$SC_UPPER $result" >> $OUTPUT
        echo "  → $SC_UPPER: $result"
    else
        echo "  $SC_UPPER: pas assez de trajectoires (n=$n_traj)"
    fi
    
    echo ""
done

echo "=== Résultats sauvegardés dans $OUTPUT ==="
echo ""
echo "Format pour publication :"
echo "  'X.XX ± Y.YY% (moyenne ± SE, n=3 trajectoires indépendantes)'"
echo ""
cat $OUTPUT
