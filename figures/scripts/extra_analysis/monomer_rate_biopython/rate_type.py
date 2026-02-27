#!/usr/bin/env python3
"""
À partir du CSV produit précédemment (orientation_{aa}.csv),
calcule les fréquences d'observation (histogramme 2D) pour deux angles :
  - theta1 : angle_z entre 2 atomes (ex. CG-CZ)
  - theta2 : normale au plan (normal_angle_z)
en fonction de la profondeur z (centroïde),
et génère un fichier freq_angle_{aa}.dat contenant :
  z_center  theta_center  count  angle_type

Usage :
    python freq_angle_analysis.py SCF [--bins_z N] [--bins_theta M]

Le script produit :
  - freq_angle_scf.dat
"""
import sys
import os
import argparse
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Histogramme 2D fréquence angle vs profondeur")
    p.add_argument("aa", nargs='?', default=None,
                   help="Code de la protéine (ex. SCF). Si absent, lit orientation_all.csv")
    p.add_argument("--bins_z", type=int, default=50,
                   help="Nombre de bins pour la profondeur z")
    p.add_argument("--bins_theta", type=int, default=50,
                   help="Nombre de bins pour les angles")
    return p.parse_args()


def load_data(aa):
    """Charge orientation_{aa}.csv ou orientation_all.csv"""
    fname = f"orientation_mono_{aa.lower() if aa else 'all'}.csv"
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"Fichier '{fname}' introuvable.")
    return pd.read_csv(fname)


def extract_records(df):
    """
    Reconstitue une table plate avec colonnes : depth_z, theta, angle_type
    - theta2 pour angle_z (entre deux atomes) et theta1 pour normal_angle_z
    On détecte automatiquement les colonnes et on associe chaque angle à sa profondeur
      (même suffixe '_trajN').
    """
    cols = df.columns.tolist()
    # colonnes angle1 (entre atomes), excluant la normale
    angle1_cols = [c for c in cols if c.endswith('_angle_z') or '_angle_z_traj' in c]
    angle1_cols = [c for c in angle1_cols if not c.startswith('normal_')]
    # colonnes angle2 (normale)
    angle2_cols = [c for c in cols if c.startswith('normal_angle_z')]
    records = []
    # theta2
    for c in angle1_cols:
        # suffixe _trajN ou vide
        suffix = c.split('angle_z')[-1]
        depth_c = 'depth_z' + suffix
        if depth_c not in df.columns:
            continue
        for z, theta in zip(df[depth_c], df[c]):
            records.append({'depth_z': z,
                            'theta': theta,
                            'angle_type': 'theta2'})
    # theta1
    for c in angle2_cols:
        suffix = c.split('normal_angle_z')[-1]
        depth_c = 'depth_z' + suffix
        if depth_c not in df.columns:
            continue
        for z, theta in zip(df[depth_c], df[c]):
            records.append({'depth_z': z,
                            'theta': theta,
                            'angle_type': 'theta1'})
    return pd.DataFrame(records)


def compute_histogram(df, bins_z, bins_theta):
    """Calcule histogramme 2D pour chaque angle_type"""
    results = []
    for angle_type in df['angle_type'].unique():
        sub = df[df['angle_type'] == angle_type]
        hist, z_edges, th_edges = np.histogram2d(
            sub['depth_z'], sub['theta'],
            bins=[bins_z, bins_theta]
        )
        # centres de bins
        z_centers  = (z_edges[:-1] + z_edges[1:]) / 2
        th_centers = (th_edges[:-1] + th_edges[1:]) / 2
        for i, zc in enumerate(z_centers):
            for j, tc in enumerate(th_centers):
                count = int(hist[i, j])
                results.append({'z_center': zc,
                                'theta_center': tc,
                                'count': count,
                                'angle_type': angle_type})
    return pd.DataFrame(results)


def main():
    args = parse_args()
    df = load_data(args.aa)
    rec = extract_records(df)
    freq = compute_histogram(rec, args.bins_z, args.bins_theta)
    out = f"freq_angle_mono_{args.aa.lower() if args.aa else 'all'}.dat"
    freq.to_csv(out, sep='\t', index=False, float_format='%.6f')
    print(f"Fréquences enregistrées dans '{out}' ({len(freq)} lignes)")

if __name__ == '__main__':
    main()







