#!/usr/bin/env python3
"""
Script pour analyser l'orientation des anneaux aromatiques à partir de fichiers de coordonnées,
avec fusion des fichiers *_coor4.dat et *_coor5.dat pour plusieurs trajectoires,
et sortie d'un fichier CSV combinant les résultats par trajectoire en colonnes distinctes.

Fonctionnalités :
- Lecture des fichiers matching un pattern donné pour chaque trajectoire.
- Fusion des frames : _coor4.dat + _coor5.dat (recalage des frames).
- Calcul du vecteur unitaire entre deux atomes spécifiés.
- Calcul du vecteur normal au plan d'un anneau aromatique (via SVD).
- Mesure de l'angle entre ces vecteurs et l'axe Z.
- Calcul de la profondeur Z (centroid du cycle aromatique).
- Sélection dynamique des trajectoires selon l'argument `aa` (SCW ou autre).
- Génération d'un DataFrame par trajectoire, renommage des colonnes avec suffixe _trajN, puis merge.
- Sauvegarde des résultats dans un fichier CSV.
"""
import os
import glob
import re
import sys
import numpy as np
import pandas as pd
from functools import reduce


def parse_coor_files(path_pattern):
    coords = {}
    max_frame = -1
    files = sorted(glob.glob(path_pattern))
    files1 = [fp for fp in files if re.search(r'_coor_mono_1\.dat$', fp)]
    files2 = [fp for fp in files if re.search(r'_coor_mono_2\.dat$', fp)]
    files3 = [fp for fp in files if re.search(r'_coor_mono_3\.dat$', fp)]
    # partie 1
    for fp in files1:
        atom = os.path.basename(fp).split('_')[0]
        with open(fp) as f:
            for line in f:
                if line.startswith('#'): continue
                frame, _, idx, resid, x, y, z = line.split()
                frame, resid = int(frame), int(resid)
                coords.setdefault(frame, {}).setdefault(resid, {})[atom] = np.array([float(x), float(y), float(z)])
                max_frame = max(max_frame, frame)
    offset = max_frame + 1
    # partie 2
    for fp in files2:
        atom = os.path.basename(fp).split('_')[0]
        with open(fp) as f:
            for line in f:
                if line.startswith('#'): continue
                frame, _, idx, resid, x, y, z = line.split()
                frame, resid = int(frame) + offset, int(resid)
                coords.setdefault(frame, {}).setdefault(resid, {})[atom] = np.array([float(x), float(y), float(z)])
                max_frame = max(max_frame, frame)
    offset = max_frame + 1
    # partie 3
    for fp in files3:
        atom = os.path.basename(fp).split('_')[0]
        with open(fp) as f:
            for line in f:
                if line.startswith('#'): continue
                frame, _, idx, resid, x, y, z = line.split()
                frame, resid = int(frame) + offset, int(resid)
                coords.setdefault(frame, {}).setdefault(resid, {})[atom] = np.array([float(x), float(y), float(z)])
    return coords


def unit_vector(p1, p2):
    v = p2 - p1
    norm = np.linalg.norm(v)
    return v / norm if norm else np.zeros_like(v)


def ring_normal(points):
    pts = np.array(points)
    centroid = pts.mean(axis=0)
    _, _, vh = np.linalg.svd(pts - centroid)
    normal = vh[-1] / np.linalg.norm(vh[-1])
    return centroid, normal


def angle_between(v, axis):
    # angle en degrés entre v et axis
    norm_v = np.linalg.norm(v)
    norm_axis = np.linalg.norm(axis)
    if norm_v == 0 or norm_axis == 0:
        return np.nan
    cosθ = np.dot(v, axis) / (norm_v * norm_axis)
    cosθ = np.clip(cosθ, -1.0, 1.0)
    return np.degrees(np.arccos(cosθ))


def analyze_coords(coords, atom_pair, ring_atoms):
    rows = []
    for frame in sorted(coords):
        for resid, atoms in coords[frame].items():
            # extraire points du cycle
            pts = [atoms[a] for a in ring_atoms if a in atoms]
            if len(pts) < 3:
                continue
            centroid, normal = ring_normal(pts)
            depth_z = centroid[2]
            # choisir l'axe de référence selon le signe de z
            axis = np.array([0,0,1]) if depth_z >= 0 else np.array([0,0,-1])

            # calcul des angles
            # angle1 entre atom_pair
            a1, a2 = atom_pair
            if a1 in atoms and a2 in atoms:
                u = unit_vector(atoms[a1], atoms[a2])
                theta2 = angle_between(u, axis)
            else:
                theta2 = np.nan
            # angle2 : normale au plan
            theta1 = angle_between(normal, axis)

            rows.append({
                'frame': frame,
                'index': resid,
                f'{a1}-{a2}_angle_z': theta2,
                'normal_angle_z': theta1,
                'depth_z': depth_z
            })
    return pd.DataFrame(rows)


def main(aa=None):
    # paramètres d'analyse
    if aa in ['SCY', 'SCYM', 'SCF']:
        atom_pair = ('CG', 'CZ')
        ring_atoms = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
      
    elif aa == 'SCW':
        atom_pair = ('CZ3', 'CE2')
        ring_atoms = ['CD2', 'CZ2', 'CZ3', 'CD1', 'CE2', 'CH2', 'CE3']
      
    elif aa in ['SCHE', 'SCHD', 'SCHP']:
        atom_pair = ('CG', 'CE1')
        ring_atoms = ['CE1', 'CD2', 'CG', 'ND1', 'NE2']
    else:
        print(f"Erreur: acide aminé '{aa}' non reconnu. Valeurs acceptées: SCY, SCYM, SCF, SCW, SCHE, SCHD, SCHP")
        return

    traj_list = [1, 2, 3]

    df_list = []
    for traj in traj_list:
        if aa:
            pattern = f"../../data/aromatics_orientation/raw_data/{aa}/traj{traj}/*_coor_mono_*.dat"
        else:
            pattern = f"../../data/aromatics_orientation/raw_data/*/traj{traj}/*_coor_mono_*.dat"
        coords = parse_coor_files(pattern)
        if not coords:
            print(f"Aucun fichier trouvé pour traj{traj} (pattern: {pattern}), je passe cette trajectoire.")
            continue
        df_traj = analyze_coords(coords, atom_pair, ring_atoms)
        if df_traj.empty:
            print(f"Aucune donnée extraite pour traj{traj}, je passe cette trajectoire.")
            continue
        # Ajouter une colonne trajectory pour identifier l'origine
        df_traj['trajectory'] = traj
        df_list.append(df_traj)

    # Concaténer verticalement toutes les trajectoires
    if not df_list:
        print("Aucun résultat à sauvegarder, fin du script.")
        return
    
    df_merged = pd.concat(df_list, ignore_index=True)
    df_merged = df_merged.sort_values(['trajectory', 'frame', 'index']).reset_index(drop=True)

    # sauvegarde
    out = f'orientation_mono_{aa or "all"}.csv'.lower()
    df_merged.to_csv(out, index=False)
    print(f"Fichier de sortie: {out}")

if __name__ == '__main__':
    aa = sys.argv[1] if len(sys.argv) > 1 else None
    main(aa)






