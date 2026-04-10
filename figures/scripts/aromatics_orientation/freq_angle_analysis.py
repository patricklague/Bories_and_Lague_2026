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

data_dir="../../data/aromatics_orientation/raw_data/total/"
out_dir="../../data/aromatics_orientation/"

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
    fname = os.path.join(data_dir, f"orientation_{aa.lower() if aa else 'all'}.csv")
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"Fichier '{fname}' introuvable.")
    return pd.read_csv(fname)


def extract_records(df, traj_suffix=None, frame_start=None, frame_end=None):
    """
    Reconstitue une table plate avec colonnes : depth_z, theta, angle_type
    - theta2 pour angle_z (entre deux atomes) et theta1 pour normal_angle_z
    Si traj_suffix est donné (ex: '_traj1'), ne prend que ces colonnes.
    Si frame_start/frame_end sont donnés, filtre les frames.
    """
    if frame_start is not None and frame_end is not None:
        df = df[(df['frame'] >= frame_start) & (df['frame'] <= frame_end)]

    cols = df.columns.tolist()

    if traj_suffix:
        angle1_cols = [c for c in cols if c.endswith('_angle_z' + traj_suffix)]
        angle1_cols = [c for c in angle1_cols if not c.startswith('normal_')]
        angle2_cols = [c for c in cols if c == 'normal_angle_z' + traj_suffix]
    else:
        angle1_cols = [c for c in cols if c.endswith('_angle_z') or '_angle_z_traj' in c]
        angle1_cols = [c for c in angle1_cols if not c.startswith('normal_')]
        angle2_cols = [c for c in cols if c.startswith('normal_angle_z')]

    records = []
    # theta2
    for c in angle1_cols:
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


def compute_histogram(df, bins_z, bins_theta, z_range=None, th_range=None):
    """Calcule histogramme 2D pour chaque angle_type avec des bins fixes"""
    results = []
    for angle_type in ['theta1', 'theta2']:
        sub = df[df['angle_type'] == angle_type]
        if len(sub) == 0:
            continue
        hist, z_edges, th_edges = np.histogram2d(
            sub['depth_z'], sub['theta'],
            bins=[bins_z, bins_theta],
            range=[z_range, th_range] if z_range and th_range else None
        )
        z_centers  = (z_edges[:-1] + z_edges[1:]) / 2
        th_centers = (th_edges[:-1] + th_edges[1:]) / 2
        results.append((angle_type, hist, z_centers, th_centers, z_edges, th_edges))
    return results


def compute_batch_mean_histogram(df_full, bins_z, bins_theta):
    """
    Calcule la moyenne des histogrammes sur 9 batches :
    3 trajectoires × 3 batches de 20000 frames (0–19999, 20000–39999, 40000–59999).
    """
    batches = [(0, 19999), (20000, 39999), (40000, 59999)]
    traj_suffixes = ['_traj1', '_traj2', '_traj3']

    # Première passe : déterminer les ranges globales pour des bins cohérents
    rec_all = extract_records(df_full)
    z_range_global = {}
    th_range_global = {}
    for angle_type in ['theta1', 'theta2']:
        sub = rec_all[rec_all['angle_type'] == angle_type]
        if len(sub) == 0:
            continue
        z_range_global[angle_type] = (sub['depth_z'].min(), sub['depth_z'].max())
        th_range_global[angle_type] = (sub['theta'].min(), sub['theta'].max())

    # Accumuler les histogrammes par batch
    hist_accum = {}  # angle_type -> list of hist arrays
    centers_cache = {}

    for traj in traj_suffixes:
        for f_start, f_end in batches:
            rec = extract_records(df_full, traj_suffix=traj, frame_start=f_start, frame_end=f_end)
            if len(rec) == 0:
                continue
            for angle_type in ['theta1', 'theta2']:
                sub = rec[rec['angle_type'] == angle_type]
                if len(sub) == 0:
                    continue
                z_r = z_range_global.get(angle_type)
                th_r = th_range_global.get(angle_type)
                hist, z_edges, th_edges = np.histogram2d(
                    sub['depth_z'], sub['theta'],
                    bins=[bins_z, bins_theta],
                    range=[z_r, th_r] if z_r and th_r else None
                )
                if angle_type not in hist_accum:
                    hist_accum[angle_type] = []
                    z_centers = (z_edges[:-1] + z_edges[1:]) / 2
                    th_centers = (th_edges[:-1] + th_edges[1:]) / 2
                    centers_cache[angle_type] = (z_centers, th_centers)
                hist_accum[angle_type].append(hist)

    # Moyenne des histogrammes + normalisation par angle_type
    results = []
    for angle_type, hists in hist_accum.items():
        mean_hist = np.mean(hists, axis=0)
        total = mean_hist.sum()
        if total > 0:
            mean_hist = mean_hist / total
        z_centers, th_centers = centers_cache[angle_type]
        for i, zc in enumerate(z_centers):
            for j, tc in enumerate(th_centers):
                count = mean_hist[i, j]
                results.append({'z_center': zc,
                                'theta_center': tc,
                                'count': count,
                                'angle_type': angle_type})
    return pd.DataFrame(results)


def main():
    args = parse_args()
    df = load_data(args.aa)
    freq = compute_batch_mean_histogram(df, args.bins_z, args.bins_theta)
    out = os.path.join(out_dir, f"freq_angle_{args.aa.lower() if args.aa else 'all'}.dat")
    freq.to_csv(out, sep='\t', index=False, float_format='%.10f')
    print(f"Fréquences enregistrées dans '{out}' ({len(freq)} lignes)")

if __name__ == '__main__':
    main()







