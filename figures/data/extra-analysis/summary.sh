#!/bin/bash
# ===============================================================
# Fusionner les 3 fichiers : moyenne + erreur standard
# ===============================================================
echo "=== Fusion des fichiers ==="

paste <(grep -v '^#' percent_traj1.dat) <(grep -v '^#' percent_traj2.dat) <(grep -v '^#' percent_traj3.dat) | \
awk 'BEGIN {
    print "#frame monomer_traj1 monomer_traj2 monomer_traj3 monomer_moy standard_error popc_alone_traj1 popc_alone_traj2 popc_alone_traj3 popc_alone_moy standard_error"
}
{
    frame=$1
    m1=$2; p1=$3
    m2=$5; p2=$6
    m3=$8; p3=$9

    # Moyennes
    moy_m=(m1+m2+m3)/3.0
    moy_p=(p1+p2+p3)/3.0

    # Erreurs standards (σ/√n)
    se_m=sqrt(((m1-moy_m)^2 + (m2-moy_m)^2 + (m3-moy_m)^2)/2.0)/sqrt(3)
    se_p=sqrt(((p1-moy_p)^2 + (p2-moy_p)^2 + (p3-moy_p)^2)/2.0)/sqrt(3)

    printf "%d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", frame, m1, m2, m3, moy_m, se_m, p1, p2, p3, moy_p, se_p
}' > percent_summary.dat

echo "=== Résumé enregistré dans percent_summary.dat ==="


