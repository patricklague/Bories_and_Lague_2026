#!/usr/bin/env bash
#
# Remplace le bloc GLPA (ou autre TAG) d’un PDB buggé par celui d’un PDB de référence
# puis ajoute une ligne TER juste après.
# Usage : ./fix_glpa.sh bugge.pdb reference.pdb > corrige.pdb
#         (option)     ./fix_glpa.sh bugge.pdb reference.pdb TAG > corrige.pdb
# bash ./fix_glpa.sh PMm-SCE/charmm-gui/namd/step5_input.pdb step5_input_PMm.pdb > PMm-SCE/charmm-gui/namd/step5_input.pdb

set -euo pipefail

bugged="$1"
reference="$2"
tag="${3:-GLPA}"        # TAG à corriger (par défaut GLPA)

awk -v ref="$reference" -v TAG="$tag" '
###############################################################################
# 1) Lecture du fichier de référence, on stocke toutes les lignes TAG souhaitées
###############################################################################
BEGIN {
    tag_col_start = 73
    tag_length    = 4
    good_count    = 0
    while ((getline line < ref) > 0) {
        rec_type = substr(line, 1, 6)
        segment  = substr(line, tag_col_start, tag_length)
        if ((rec_type == "ATOM  " || rec_type == "HETATM") && segment == TAG) {
            good[++good_count] = line          # lignes GLPA correctes
        }
    }
    close(ref)
}

###############################################################################
# 2) Traitement du fichier buggé
###############################################################################
{
    rec_type = substr($0, 1, 6)
    segment  = substr($0, tag_col_start, tag_length)

    # a) Première ligne GLPA buggée rencontrée : on injecte le bloc correct
    if ((rec_type == "ATOM  " || rec_type == "HETATM") && segment == TAG) {
        if (!printed_good) {
            for (i = 1; i <= good_count; i++) print good[i]
            print "TER"                     # NEW — ajoute systématiquement TER
            printed_good = 1
        }
        next                                # on saute la ligne buggée
    }

    # b) Tout le reste est recopié tel quel
    print
}
' "$bugged"

