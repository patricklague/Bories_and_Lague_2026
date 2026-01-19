#!/bin/bash
# ===============================================================
# Script d'analyse par lots OPTIMISÉ avec parallélisation
# Utilise GNU Parallel pour traiter plusieurs lots simultanément
# et configure VMD pour utiliser le multithreading
# ===============================================================

# Configuration de la parallélisation
NUM_PARALLEL_JOBS=4      # Nombre de lots à traiter en parallèle (ajuster selon votre RAM)
VMD_NUM_THREADS=4        # Nombre de threads pour VMD

# Variables d'environnement pour VMD
export VMDNUMCPUS=$VMD_NUM_THREADS
export VMDFORCECPUCOUNT=$VMD_NUM_THREADS

# Vérifier si GNU Parallel est installé
if ! command -v parallel &> /dev/null; then
    echo "GNU Parallel n'est pas installé. Installation recommandée:"
    echo "  macOS: brew install parallel"
    echo "  Ubuntu: sudo apt install parallel"
    echo "Passage en mode séquentiel..."
    USE_PARALLEL=false
else
    USE_PARALLEL=true
    # Désactiver le message de citation de GNU Parallel
    mkdir -p ~/.parallel
    touch ~/.parallel/will-cite
fi

# Liste des analogues
aafile=("SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCQ" "SCN" "SCF" "SCW" "GLYD" "SCP" "SCCM" "SCYM" "SCE" "SCEN" "SCD" "SCDN" "SCK" "SCKN" "SCR" "SCRN" "SCHE" "SCHD" "SCHP")
#aafile=("SCI")  # Pour tester avec un seul analogue

# Paramètres des lots
BATCH_SIZE=50       # Nombre de sections par lot
SECTION_START=400   # Première section
SECTION_END=1000    # Dernière section

# Répertoire de travail temporaire
WORK_DIR=$(pwd)
TEMP_BASE="/tmp/vmd_analysis_$$"
mkdir -p "$TEMP_BASE"

# Fonction pour traiter un lot individuel
process_batch() {
    local aa=$1
    local t=$2
    local start=$3
    local end=$4
    local batch_num=$5
    local DIR=$6
    local PSFFILE=$7
    local PDBFILE=$8
    
    # Créer un répertoire temporaire unique pour ce lot
    local TEMP_DIR="$TEMP_BASE/${aa}_traj${t}_batch${batch_num}"
    mkdir -p "$TEMP_DIR"
    cd "$TEMP_DIR"
    
    # Copier les scripts nécessaires
    cp "$WORK_DIR/indexNoWat.vmd" .
    cp "$WORK_DIR/center.vmd" .
    cp "$WORK_DIR/indexNoMemb.vmd" .
    cp "$WORK_DIR/densityProfiles-aa.vmd" .
    cp "$WORK_DIR/average_profiles.py" .
    cp "$WORK_DIR/input.psf" .
    
    # Construire la liste des fichiers DCD pour ce lot
    local files=""
    for (( i=$start; i<=$end; i++ )); do
        files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    
    # Retirer l'eau
    cat ./indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp1.vmd
    vmd -dispdev text -nt $VMD_NUM_THREADS -e temp1.vmd >& /dev/null
    
    local PSFFILE2="./analog-$aa.psf"
    local PDBFILE2="./analog-$aa.pdb"
    
    # Concaténer les DCD avec stride=1
    catdcd -o trajectory.dcd -stride 1 -i findexfile.ind $files
    
    # Centrer la trajectoire
    local DCDFILE="trajectory.dcd"
    cat ./center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE=$DCDFILE= | sed s=AAA=$aa= > temp2.vmd
    vmd -dispdev text -nt $VMD_NUM_THREADS -e ./temp2.vmd >& /dev/null
    rm trajectory.dcd
    
    # Supprimer la membrane
    cat ./indexNoMemb.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=PDBFILE=$PDBFILE2= | sed s=AAA=$aa= > temp3.vmd
    vmd -dispdev text -nt $VMD_NUM_THREADS -e ./temp3.vmd >& /dev/null
    
    local PSFFILE3="./analog-$aa.psf"
    catdcd -o trajectory.dcd -stride 1 -i findexfile.ind centered.dcd
    
    # Préparer et exécuter le script TCL pour les profils de densité
    cat ./densityProfiles-aa.vmd | sed s=AAA=$aa= | sed s=PSFFILE=$PSFFILE3= | sed s=DCDFILE="./trajectory.dcd"= > solutes.vmd
    vmd -dispdev text -nt $VMD_NUM_THREADS -e solutes.vmd >& /dev/null
    
    # Calculer les moyennes
    python ./average_profiles.py dens_mono_frame_*.dat > dens_mono_traj_${t}_batch_${batch_num}.dat
    python ./average_profiles.py dens_nonmono_frame_*.dat > dens_nonmono_traj_${t}_batch_${batch_num}.dat
    
    # Copier les résultats vers le répertoire de travail
    cp dens_mono_traj_${t}_batch_${batch_num}.dat "$WORK_DIR/"
    cp dens_nonmono_traj_${t}_batch_${batch_num}.dat "$WORK_DIR/"
    
    # Nettoyer
    cd "$WORK_DIR"
    rm -rf "$TEMP_DIR"
    
    echo "  Batch $batch_num completed for $aa traj $t"
}

# Exporter la fonction pour GNU Parallel
export -f process_batch
export TEMP_BASE WORK_DIR VMD_NUM_THREADS

# Boucle principale
for aa in "${aafile[@]}"; do
    # Déterminer le répertoire source
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
    traj=(1 2 3)
    
    if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
        DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
    elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
        traj=(4 2 3)
    fi
    
    PSFFILE="$DIR/charmm-gui/namd/step5_input.psf"
    PDBFILE="$DIR/charmm-gui/namd/step5_input.pdb"
    
    # Copier le fichier PSF
    cp ../../../supp_files/input_systems/${aa,,}/step5_input.psf ./input.psf
    
    for t in "${traj[@]}"; do
        echo "=============================================="
        echo "Processing $aa : trajectory $t"
        echo "=============================================="
        
        DATA="$DIR/analyses/traj$t/data"
        
        if [ "$USE_PARALLEL" = true ]; then
            # Mode parallèle avec GNU Parallel
            # Créer un fichier temporaire avec la liste des lots
            batch_file=$(mktemp)
            batch_num=0
            for (( start=$SECTION_START; start<=$SECTION_END; start+=$BATCH_SIZE )); do
                end=$((start + BATCH_SIZE - 1))
                if [ $end -gt $SECTION_END ]; then
                    end=$SECTION_END
                fi
                echo "$aa $t $start $end $batch_num $DIR $PSFFILE $PDBFILE" >> "$batch_file"
                batch_num=$((batch_num + 1))
            done
            
            # Exécuter en parallèle
            cat "$batch_file" | parallel -j $NUM_PARALLEL_JOBS --colsep ' ' \
                process_batch '{1}' '{2}' '{3}' '{4}' '{5}' '{6}' '{7}' '{8}'
            
            # Nettoyer le fichier temporaire
            rm -f "$batch_file"
        else
            # Mode séquentiel (fallback)
            batch_num=0
            for (( start=$SECTION_START; start<=$SECTION_END; start+=$BATCH_SIZE )); do
                end=$((start + BATCH_SIZE - 1))
                if [ $end -gt $SECTION_END ]; then
                    end=$SECTION_END
                fi
                
                echo "  Processing batch $batch_num: sections $start to $end"
                process_batch "$aa" "$t" "$start" "$end" "$batch_num" "$DIR" "$PSFFILE" "$PDBFILE"
                
                batch_num=$((batch_num + 1))
            done
        fi
        
        # Combiner tous les résultats des lots
        python ./average_profiles.py dens_mono_traj_${t}_batch*.dat > dens_mono_avg.dat
        python ./average_profiles.py dens_nonmono_traj_${t}_batch*.dat > dens_nonmono_avg.dat
        ./combine_profiles.sh dens_mono_avg.dat dens_nonmono_avg.dat > dens_all_avg.dat
        
        mv -f dens_all_avg.dat $DATA/densityProfiles/profile-aa-600-mono-2.dat
        mkdir -p contacts/${aa,,}
        mv dens_*mono_traj_${t}_batch*.dat contacts/${aa,,}/
    done
done

# Nettoyer le répertoire temporaire
rm -rf "$TEMP_BASE"

echo "=============================================="
echo "All analyses completed!"
echo "=============================================="
