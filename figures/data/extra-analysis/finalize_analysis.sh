#!/bin/bash
# ===============================================================
# Script pour finaliser l'analyse après tous les lots
# Fusionne les fichiers de lots et calcule les statistiques finales
#
# Usage: ./finalize_analysis.sh <ANALOGUE>
# ===============================================================

if [ $# -lt 1 ]; then
    echo "Usage: $0 <ANALOGUE>"
    exit 1
fi

aa=${1^^}  # Majuscules
aa_lower=${aa,,}  # Minuscules

BATCH_DIR="contacts/${aa_lower}/batches"
FRAMES_PER_SECTION=100
BATCH_SIZE=10

echo "=============================================="
echo "Finalizing analysis for $aa"
echo "=============================================="

# Vérifier que le répertoire existe
if [ ! -d "$BATCH_DIR" ]; then
    echo "Error: Directory $BATCH_DIR not found!"
    exit 1
fi

# --- Fusionner les fichiers de pourcentage pour chaque trajectoire ---
for traj_num in 1 2 3
do
    echo ""
    echo "Processing trajectory $traj_num..."
    
    OUTPUT_FILE="contacts/${aa_lower}/percent_traj${traj_num}.dat"
    echo "#frame monomers_% popc_isolated_%" > "$OUTPUT_FILE"
    
    frame_offset=0
    batch_num=0
    
    while true
    do
        BATCH_FILE="${BATCH_DIR}/percent_traj${traj_num}_batch${batch_num}.dat"
        
        if [ ! -f "$BATCH_FILE" ]; then
            break
        fi
        
        echo "  Merging batch $batch_num..."
        
        # Ajouter les frames avec le bon offset
        awk -v offset=$frame_offset 'NR>1 && !/^#/ && NF>=3 {
            printf "%d %.6f %.6f\n", ($1 + offset), $2, $3
        }' "$BATCH_FILE" >> "$OUTPUT_FILE"
        
        # Compter les frames
        frames_in_batch=$(grep -v '^#' "$BATCH_FILE" | grep -c '.')
        frame_offset=$((frame_offset + frames_in_batch))
        
        batch_num=$((batch_num + 1))
    done
    
    echo "  -> Total frames: $frame_offset"
    echo "  -> Output: $OUTPUT_FILE"
done

# --- Fusionner les fichiers de persistance ---
for traj_num in 1 2 3
do
    echo ""
    echo "Merging persistence data for trajectory $traj_num..."
    
    OUTPUT_PERS="contacts/${aa_lower}/persistence_traj${traj_num}.dat"
    
    # Pour la persistance, on doit agréger les statistiques par resid
    # On collecte toutes les durées de tous les lots et on recalcule
    
    # Créer un fichier temporaire pour collecter toutes les données
    TMP_FILE=$(mktemp)
    
    batch_num=0
    while true
    do
        PERS_FILE="${BATCH_DIR}/persistence_traj${traj_num}_batch${batch_num}.dat"
        
        if [ ! -f "$PERS_FILE" ]; then
            break
        fi
        
        # Extraire les données (ignorer les commentaires)
        grep -v '^#' "$PERS_FILE" >> "$TMP_FILE"
        
        batch_num=$((batch_num + 1))
    done
    
    # Agréger par resid: pour chaque résidu, calculer la moyenne pondérée
    echo "#resid average_duration max_duration persistent_contacts" > "$OUTPUT_PERS"
    
    # Trier par resid et agréger
    sort -n "$TMP_FILE" | awk '
    BEGIN { prev_resid = -1 }
    {
        resid = $1
        avg_dur = $2
        max_dur = $3
        contacts = $4
        
        if (resid != prev_resid && prev_resid != -1) {
            # Écrire le résidu précédent
            if (total_contacts > 0) {
                final_avg = total_weighted_dur / total_contacts
            } else {
                final_avg = 0
            }
            printf "%d %.3f %d %d\n", prev_resid, final_avg, global_max, total_contacts
            
            # Reset
            total_weighted_dur = 0
            total_contacts = 0
            global_max = 0
        }
        
        # Accumuler
        total_weighted_dur += avg_dur * contacts
        total_contacts += contacts
        if (max_dur > global_max) global_max = max_dur
        
        prev_resid = resid
    }
    END {
        if (prev_resid != -1) {
            if (total_contacts > 0) {
                final_avg = total_weighted_dur / total_contacts
            } else {
                final_avg = 0
            }
            printf "%d %.3f %d %d\n", prev_resid, final_avg, global_max, total_contacts
        }
    }' >> "$OUTPUT_PERS"
    
    rm -f "$TMP_FILE"
    echo "  -> Output: $OUTPUT_PERS"
done

# --- Calculer les statistiques finales (moyenne des 3 trajectoires) ---
echo ""
echo "Computing final statistics..."

cd "contacts/${aa_lower}"

# Vérifier que tous les fichiers existent
if [ -f "percent_traj1.dat" ] && [ -f "percent_traj2.dat" ] && [ -f "percent_traj3.dat" ]; then
    bash ../../summary_percent_full.sh
    mv percent_summary.dat "percent_${aa_lower}.dat"
    echo "  -> percent_${aa_lower}.dat created"
else
    echo "  Warning: Not all percent_trajX.dat files found"
fi

if [ -f "persistence_traj1.dat" ] && [ -f "persistence_traj2.dat" ] && [ -f "persistence_traj3.dat" ]; then
    bash ../../summary_persistence_full.sh
    mv persistence_summary.dat "persistence_${aa_lower}.dat"
    echo "  -> persistence_${aa_lower}.dat created"
else
    echo "  Warning: Not all persistence_trajX.dat files found"
fi

cd ../..

echo ""
echo "=============================================="
echo "Finalization complete for $aa"
echo "Results in: contacts/${aa_lower}/"
echo "=============================================="
