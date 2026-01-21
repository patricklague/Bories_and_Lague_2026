#!/bin/bash
# ===============================================================
# Script pour fusionner les fichiers de lots d'un même trajectoire
# Usage: merge_batches.sh <analogue_lowercase> <traj_num> <num_batches>
# ===============================================================

AA=$1          # ex: "sci"
TRAJ_NUM=$2    # ex: 1, 2 ou 3
NUM_BATCHES=$3 # Nombre total de lots

BATCH_DIR="contacts/${AA}/batches"
OUTPUT_FILE="contacts/${AA}/percent_traj${TRAJ_NUM}.dat"
OUTPUT_PERS="contacts/${AA}/persistence_traj${TRAJ_NUM}.dat"

echo "Merging $NUM_BATCHES batches for $AA trajectory $TRAJ_NUM..."

# --- Fusion des fichiers de pourcentage ---
# Chaque lot a 1000 frames (10 sections * 100 frames/section)
# On doit recalculer les numéros de frame globaux

echo "#frame monomers_% popc_isolated_%" > $OUTPUT_FILE

frame_offset=0
for (( b=0; b<$NUM_BATCHES; b++ ))
do
  BATCH_FILE="${BATCH_DIR}/percent_traj${TRAJ_NUM}_batch${b}.dat"
  
  if [ -f "$BATCH_FILE" ]; then
    # Ajouter les frames avec offset, ignorer l'en-tête
    awk -v offset=$frame_offset 'NR>1 && !/^#/ {
      print ($1 + offset), $2, $3
    }' "$BATCH_FILE" >> $OUTPUT_FILE
    
    # Compter les frames dans ce lot pour l'offset suivant
    frames_in_batch=$(grep -v '^#' "$BATCH_FILE" | wc -l)
    frame_offset=$((frame_offset + frames_in_batch))
  else
    echo "Warning: $BATCH_FILE not found!"
  fi
done

echo "  Merged percent data: $frame_offset total frames -> $OUTPUT_FILE"

# --- Fusion des fichiers de persistance ---
# Pour la persistance, on doit faire attention car les contacts peuvent
# traverser les frontières de lots. Pour l'instant, on concatène simplement
# les statistiques moyennes.

echo "#batch avg_duration max_duration contacts_count" > $OUTPUT_PERS

for (( b=0; b<$NUM_BATCHES; b++ ))
do
  PERS_FILE="${BATCH_DIR}/persistence_traj${TRAJ_NUM}_batch${b}.dat"
  
  if [ -f "$PERS_FILE" ]; then
    # Extraire les statistiques globales du lot
    # Le format exact dépend de check_contact_persistency.tcl
    # Pour l'instant, on copie les statistiques de chaque lot
    tail -1 "$PERS_FILE" | awk -v b=$b '{print b, $0}' >> $OUTPUT_PERS
  fi
done

echo "  Merged persistence data -> $OUTPUT_PERS"
echo "Done merging trajectory $TRAJ_NUM for $AA"
