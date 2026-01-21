#!/bin/bash
# ===============================================================
# Script d'analyse par lots - Version avec contrôle fin
# Permet de traiter un analogue spécifique ou de reprendre une analyse
# 
# Usage:
#   ./launch_analysis_batched_single.sh <ANALOGUE> [START_BATCH] [END_BATCH]
#
# Exemples:
#   ./launch_analysis_batched_single.sh SCI           # Tous les lots
#   ./launch_analysis_batched_single.sh SCI 0 9       # Lots 0-9 seulement
#   ./launch_analysis_batched_single.sh SCI 10 19     # Lots 10-19 seulement
# ===============================================================

if [ $# -lt 1 ]; then
    echo "Usage: $0 <ANALOGUE> [START_BATCH] [END_BATCH]"
    echo "  ANALOGUE: nom de l'analogue (ex: SCI, SCY, GLYD)"
    echo "  START_BATCH: premier lot à traiter (défaut: 0)"
    echo "  END_BATCH: dernier lot à traiter (défaut: tous)"
    exit 1
fi

aa=${1^^}  # Convertir en majuscules
START_BATCH=${2:-0}
END_BATCH=${3:--1}  # -1 = jusqu'à la fin

# Paramètres des lots
BATCH_SIZE=10       # 10 sections par lot = 1000 frames
SECTION_START=400
SECTION_END=1000
FRAMES_PER_SECTION=100  # Chaque fichier DCD contient 100 frames

# Calculer le nombre total de lots
TOTAL_SECTIONS=$((SECTION_END - SECTION_START + 1))  # 601 sections
TOTAL_BATCHES=$(( (TOTAL_SECTIONS + BATCH_SIZE - 1) / BATCH_SIZE ))  # ~61 lots

if [ $END_BATCH -eq -1 ]; then
    END_BATCH=$((TOTAL_BATCHES - 1))
fi

echo "=============================================="
echo "Analogue: $aa"
echo "Batch size: $BATCH_SIZE sections ($((BATCH_SIZE * FRAMES_PER_SECTION)) frames)"
echo "Total batches: $TOTAL_BATCHES"
echo "Processing batches: $START_BATCH to $END_BATCH"
echo "=============================================="

# Créer les répertoires
mkdir -p contacts/${aa,,}
mkdir -p contacts/${aa,,}/batches

# Déterminer le répertoire source
DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
traj=(1 2 3)
if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
    traj=(4 2 3)
fi

echo "Source directory: $DIR"
echo "Trajectories: ${traj[@]}"

# Copier le fichier PSF
cp ../../../supp_files/input_systems/${aa,,}/step5_input.psf ./input.psf

for t in "${traj[@]}"
do
    echo ""
    echo "=============================================="
    echo "Processing trajectory $t"
    echo "=============================================="
    
    # Numéro de trajectoire pour les noms de fichiers
    if [ $t = 4 ]; then
        traj_num=1
    else
        traj_num=$t
    fi
    
    # Boucle sur les lots demandés
    for (( batch_num=$START_BATCH; batch_num<=$END_BATCH; batch_num++ ))
    do
        # Calculer les sections pour ce lot
        start_section=$((SECTION_START + batch_num * BATCH_SIZE))
        end_section=$((start_section + BATCH_SIZE - 1))
        
        if [ $start_section -gt $SECTION_END ]; then
            echo "Batch $batch_num: no more sections to process"
            break
        fi
        
        if [ $end_section -gt $SECTION_END ]; then
            end_section=$SECTION_END
        fi
        
        echo "  Batch $batch_num: sections $start_section to $end_section"
        
        # Vérifier si le fichier existe déjà (pour reprendre)
        OUT_FILE="contacts/${aa,,}/batches/percent_traj${traj_num}_batch${batch_num}.dat"
        if [ -f "$OUT_FILE" ]; then
            echo "    -> Already processed, skipping"
            continue
        fi
        
        # Construire la liste des fichiers DCD
        files=""
        for (( i=$start_section; i<=$end_section; i++ ))
        do
            files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
        done
        
        # Concaténer avec stride=1
        echo "    Creating trajectory.dcd..."
        catdcd -o trajectory.dcd -stride 1 $files
        
        # Analyse des pourcentages de monomères
        echo "    Running monomer analysis..."
        sed -r "s/AAA/$aa/g" get_purcentage.tcl > temp.tcl
        vmd -dispdev text -e temp.tcl
        mv monomers_percent.dat "$OUT_FILE"
        
        # Analyse de persistance
        echo "    Running persistence analysis..."
        sed -r "s/AAA/$aa/g" check_contact_persistency.tcl > temp.tcl
        vmd -dispdev text -e temp.tcl
        mv contact_persistence.dat "contacts/${aa,,}/batches/persistence_traj${traj_num}_batch${batch_num}.dat"
        
        # Nettoyer
        rm -f trajectory.dcd temp.tcl
        
        echo "    -> Done"
    done
done

rm -f input.psf

echo ""
echo "=============================================="
echo "Batch processing completed for $aa"
echo "Results saved in contacts/${aa,,}/batches/"
echo ""
echo "To merge batches, run:"
echo "  ./finalize_analysis.sh $aa"
echo "=============================================="
