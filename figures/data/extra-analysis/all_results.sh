#!/bin/bash
#!/bin/bash
# ===============================================================
# Combine les fichiers percent_monomer_<mol>.dat pour plusieurs molécules
# Produit un résumé global : mol_name, moyennes et erreurs standard
# ===============================================================
# Usage :
#   ./summary_all_molecules.sh mol1 mol2 mol3 ...
# (ou avec un tableau aafile[@])
# ===============================================================

# Si tu préfères définir la liste directement ici :
#aafile=("SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW")
aafile=("SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW")
# Ou bien la passer en argument :
# aafile=("$@")

# Fichier de sortie global
outfile="monomer_all.dat"
#echo "#mol_name monomer_rate_per_traj standard_error popc_alone_rate_per_traj standard_error" > $outfile
echo "#mol_name monomer_rate_per_traj standard_error" > $outfile

# Boucle sur chaque molécule
for aa in "${aafile[@]}"; do
    file=contacts/${aa,,}/"percent_${aa,,}.dat"
    if [ ! -f "$file" ]; then
        echo "⚠️  Fichier $file introuvable — ignoré"
        continue
    fi

    # On ignore la ligne de header (#)
    # On calcule la moyenne et erreur standard de toutes les frames pour cette molécule
    awk -v mol="$aa" '
        BEGIN {sum_m=0; sum_m2=0; sum_p=0; sum_p2=0; n=0}
        !/^#/ {
            m=$5;  # moyenne monomers_% (colonne 5 du summary précédent)
            se_m=$6; # erreur standard monomers
            p=$10; # moyenne popc_alone_% (colonne 10)
            se_p=$11; # erreur standard popc
            sum_m+=m; sum_m2+=m*m;
            sum_p+=p; sum_p2+=p*p;
            n++
        }
        END {
            if (n>0) {
                mean_m = sum_m/n;
                sd_m = sqrt((sum_m2 - n*mean_m^2)/(n-1));
                se_m = sd_m/sqrt(n);

                mean_p = sum_p/n;
                sd_p = sqrt((sum_p2 - n*mean_p^2)/(n-1));
                se_p = sd_p/sqrt(n);

                printf "%s %.3f %.3f \n", mol, mean_m, se_m
            }
        }' "$file" >> $outfile
done
#printf "%s %.3f %.3f %.3f %.3f\n", mol, mean_m, se_m, mean_p, se_p

echo "✅ Résumé global enregistré dans $outfile"


# Fichier de sortie
outfile="persistence_all.dat"
echo "#mol_name avg_per_traj standard_error max_max_duration nb_total_of_persistent_contacts nb_frames" > "$outfile"

# Boucle sur chaque molécule
for aa in "${aafile[@]}"; do
    file="contacts/${aa,,}/persistence_${aa,,}.dat"

    if [ ! -f "$file" ]; then
        echo "⚠️  Fichier $file introuvable — ignoré"
        continue
    fi

    # Lecture et analyse avec awk
    awk -v mol="$aa" '
        BEGIN {
            sum_avg = 0; sum_avg2 = 0; n = 0;
            max_maxdur = 0;
            total_persistent = 0;
            nb_frames = 0;
        }
        # Ignore les lignes de header
    !/^#/ {
        # Colonnes issues du summary combiné :
        # 1: resid
        # 2-4: avg_traj1 avg_traj2 avg_traj3
        # 5-7: max_traj1 max_traj2 max_traj3
        # 8-10: pers_traj1 pers_traj2 pers_traj3

        avg1 = $2; avg2 = $3; avg3 = $4;
        max1 = $5; max2 = $6; max3 = $7;
        pers1 = $8; pers2 = $9; pers3 = $10;

        # Moyenne des trois trajectoires
        mean_avg = (avg1 + avg2 + avg3) / 3.0;

        # Somme et somme des carrés pour calculer moyenne et erreur standard ensuite
        sum_avg += mean_avg;
        sum_avg2 += mean_avg * mean_avg;

        # Maximum global sur les trois trajectoires
        maxdur = max1;
        if (max2 > maxdur) maxdur = max2;
        if (max3 > maxdur) maxdur = max3;
        if (maxdur > max_maxdur) max_maxdur = maxdur;

        # Total des contacts persistants sur les trois traj
        total_persistent += pers1 + pers2 + pers3;

        n++;
    }
        /^#TOTAL/ {
            nb_frames = $5;
        }
        END {
            if (n > 1) {
                mean = sum_avg / n;
                sd = sqrt((sum_avg2 - n*mean*mean)/(n-1));
                se = sd / sqrt(n);
            } else {
                mean = (n == 1) ? sum_avg : 0;
                se = 0;
            }
            printf "%s %.3f %.3f %.1f %d %d\n", mol, mean, se, max_maxdur, total_persistent, nb_frames;
        }
    ' "$file" >> "$outfile"
done

echo "✅ Résumé global enregistré dans $outfile"



