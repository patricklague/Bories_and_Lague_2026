#!/usr/bin/env bash
set -euo pipefail          # Stoppe sur erreur, variable indéfinie ou pipe cassé
shopt -s nullglob          # Ignore la boucle si aucun fichier ne correspond

for file in *-aa/*-*/charmm-gui/namd/step7_production.inp; do
  printf '▶️  %s\n' "$file"

  # --------------------------------------------------------------------
  # 1) Sauvegarde puis remplacement de dcdfreq
  # --------------------------------------------------------------------
  [[ -e "${file}.bak" ]] || cp -- "$file" "${file}.bak"

  sed -i -E 's/^(dcdfreq[[:space:]]+)[0-9]+;/\15000;/' "$file"

  # ------------------------------------------------------------------
  # 2) détection de SCHP / SCHE / SCP / SCEN n’importe où dans le chemin
  # ------------------------------------------------------------------
  dir_path=${file%/*}   # enlève « /step7_production.inp »
  case "$dir_path" in
    *SCHP*|*SCHE*|*SCP*|*SCEN*)
      if ! grep -q '^[[:space:]]*parameters[[:space:]]\+toppar/top_all36_sidechains\.str' "$file"; then
        sed -i '/^[[:space:]]*source[[:space:]]/i parameters              toppar/top_all36_sidechains.str' "$file"
        echo "    ➕ Ligne side‑chains ajoutée"
      else
        echo "    ↩︎ Ligne présente, rien à faire"
      fi
      ;;
  esac
done

echo "✅ Modifications terminées (sauvegardes *.bak disponibles)."
