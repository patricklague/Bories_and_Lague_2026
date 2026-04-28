#!/usr/bin/env bash
#
# Script pour extraire, depuis freq_angle_{AA}.dat,
# les angles (theta_center) de plus forte fréquence
# pour z_center entre 0 et 20 Å, pour chaque angle_type.
# Usage :
#   ./extract_top_angles.sh SCF    # pour un AA particulier
#   ./extract_top_angles.sh all    # pour 'all'

AA=${1:-all}
IN="../../data/aromatics_orientation/freq_angle_${AA,,}.dat"

if [[ ! -f "$IN" ]]; then
  echo "Fichier '$IN' introuvable." >&2
  exit 1
fi

echo "# Résultats pour AA = $AA (10 Å ≤ z ≤ 20 Å)"
echo "angle_type\tz_center\ttheta_center\tcount"

awk -F '\t' 'NR>1 && $1 >= 10 && $1 <= 20 {
    # $1=z_center, $2=theta_center, $3=count, $4=angle_type
    if ($4=="theta1" && $3>max1) { max1=$3; z1=$1; t1=$2 }
    if ($4=="theta2" && $3>max2) { max2=$3; z2=$1; t2=$2 }
}
END {
    if (z1!="") printf("theta1\t%.3f\t%.3f\t%d\n", z1, t1, max1)
    if (z2!="") printf("theta2\t%.3f\t%.3f\t%d\n", z2, t2, max2)
}' "$IN"

