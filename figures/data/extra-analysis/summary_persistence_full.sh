#!/bin/bash
# ===============================================================
# Fusionner les fichiers de persistance des 3 trajectoires
# Pour les données complètes (60,000 frames par trajectoire)
# ===============================================================
echo "=== Fusion des fichiers de persistance (full data) ==="

# Vérifier que les fichiers existent
if [ ! -f "persistence_traj1.dat" ] || [ ! -f "persistence_traj2.dat" ] || [ ! -f "persistence_traj3.dat" ]; then
    echo "Error: Missing persistence_trajX.dat files!"
    exit 1
fi

# Format d'entrée: #resid average_duration max_duration persistent_contacts
# On calcule la moyenne et erreur standard pour chaque résidu

echo "#resid avg_dur_moy avg_dur_se max_dur_moy max_dur_se contacts_moy contacts_se" > persistence_summary.dat

# Extraire la liste des resids (communs aux 3 trajectoires)
# En supposant que les resids sont les mêmes dans les 3 fichiers
paste <(grep -v '^#' persistence_traj1.dat) \
      <(grep -v '^#' persistence_traj2.dat) \
      <(grep -v '^#' persistence_traj3.dat) | \
awk '{
    resid=$1
    
    # Trajectory 1
    avg1=$2; max1=$3; cnt1=$4
    # Trajectory 2
    avg2=$6; max2=$7; cnt2=$8
    # Trajectory 3
    avg3=$10; max3=$11; cnt3=$12

    # Moyennes
    moy_avg=(avg1+avg2+avg3)/3.0
    moy_max=(max1+max2+max3)/3.0
    moy_cnt=(cnt1+cnt2+cnt3)/3.0

    # Erreurs standards (σ/√n)
    se_avg=sqrt(((avg1-moy_avg)^2 + (avg2-moy_avg)^2 + (avg3-moy_avg)^2)/2.0)/sqrt(3)
    se_max=sqrt(((max1-moy_max)^2 + (max2-moy_max)^2 + (max3-moy_max)^2)/2.0)/sqrt(3)
    se_cnt=sqrt(((cnt1-moy_cnt)^2 + (cnt2-moy_cnt)^2 + (cnt3-moy_cnt)^2)/2.0)/sqrt(3)

    printf "%d %.3f %.3f %.3f %.3f %.3f %.3f\n", resid, moy_avg, se_avg, moy_max, se_max, moy_cnt, se_cnt
}' >> persistence_summary.dat

echo "=== Résumé de persistance enregistré dans persistence_summary.dat ==="
