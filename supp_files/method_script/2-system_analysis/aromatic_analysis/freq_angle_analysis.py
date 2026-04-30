#!/usr/bin/env python3
"""
From the CSV per trajectory produced by `ring_orientation_analysis.py`
(orientation_{aa}_traj{1,2,3}.csv), calculate the observation frequencies
(2D histogram) for two angles:
  - theta1: angle between the normal to the ring plane and the Z axis (normal_angle_z)
  - theta2: angle between the vector between 2 atoms (e.g., CG-CZ) and the Z axis
as a function of the depth z (ring centroid),
and generates a file freq_angle_{aa}.dat containing :
  z_center  theta_center  count  angle_type

Layout of the input data:
  figures/data/aromatics_orientation/raw_data/total/{AA}/vector_orientations/
      orientation_{aa}_traj1.csv
      orientation_{aa}_traj2.csv
      orientation_{aa}_traj3.csv

Each CSV has columns:
  frame, index, {a1}-{a2}_angle_z, normal_angle_z, depth_z

Output:
  figures/data/aromatics_orientation/freq_angle_{aa}.dat

Usage :
    python freq_angle_analysis.py SCF [--bins_z N] [--bins_theta M]
"""
import os
import argparse
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Repo layout (script lives in
# supp_files/method_script/2-system_analysis/aromatic_analysis/, so 4 ".." to
# reach the repo root).
# ---------------------------------------------------------------------------
DATA_ROOT = os.path.join(
    "..", "..", "..", "..",
    "figures", "data", "aromatics_orientation", "raw_data", "total",
)
OUT_DIR = os.path.join(
    "..", "..", "..", "..",
    "figures", "data", "aromatics_orientation"
)

# 3 trajectories × 3 sections of 20000 frames -> 9 batches
BATCHES = [(0, 19999), (20000, 39999), (40000, 59999)]
TRAJS = (1, 2, 3)


def parse_args():
    p = argparse.ArgumentParser(description="Histogramme 2D fréquence angle vs profondeur")
    p.add_argument("aa", help="Code de l'analogue (ex. SCF, SCY, SCW)")
    p.add_argument("--bins_z", type=int, default=50,
                   help="Nombre de bins pour la profondeur z")
    p.add_argument("--bins_theta", type=int, default=50,
                   help="Nombre de bins pour les angles")
    return p.parse_args()


def load_traj(aa: str, traj: int) -> pd.DataFrame:
    """Read orientation_{aa}_traj{traj}.csv for one trajectory."""
    fname = os.path.join(
        DATA_ROOT, aa, "vector_orientations",
        f"orientation_{aa.lower()}_traj{traj}.csv",
    )
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"File '{fname}' not found.")
    return pd.read_csv(fname)


def extract_records(df: pd.DataFrame,
                    frame_start: int | None = None,
                    frame_end: int | None = None) -> pd.DataFrame:
    """Reshape one per-trajectory df into long format: depth_z, theta, angle_type.

    - theta1 = normal_angle_z
    - theta2 = the *_angle_z column that is NOT 'normal_angle_z'
    """
    if frame_start is not None and frame_end is not None:
        df = df[(df["frame"] >= frame_start) & (df["frame"] <= frame_end)]

    if df.empty or "depth_z" not in df.columns:
        return pd.DataFrame(columns=["depth_z", "theta", "angle_type"])

    angle1_cols = [c for c in df.columns
                   if c.endswith("_angle_z") and c != "normal_angle_z"]

    records = []
    # theta2: pair-vector angle
    for c in angle1_cols:
        sub = df[["depth_z", c]].rename(columns={c: "theta"})
        sub["angle_type"] = "theta2"
        records.append(sub)
    # theta1: normal angle
    if "normal_angle_z" in df.columns:
        sub = df[["depth_z", "normal_angle_z"]].rename(
            columns={"normal_angle_z": "theta"})
        sub["angle_type"] = "theta1"
        records.append(sub)

    if not records:
        return pd.DataFrame(columns=["depth_z", "theta", "angle_type"])
    return pd.concat(records, ignore_index=True).dropna()


def compute_batch_mean_histogram(traj_dfs: dict, bins_z: int, bins_theta: int) -> pd.DataFrame:
    """Average 2D histograms over 3 trajectories x 3 frame batches = 9 batches."""
    # Pass 1: global ranges per angle_type for consistent bins
    rec_all = pd.concat([extract_records(df) for df in traj_dfs.values()],
                        ignore_index=True)
    z_range_global = {}
    th_range_global = {}
    for angle_type in ("theta1", "theta2"):
        sub = rec_all[rec_all["angle_type"] == angle_type]
        if sub.empty:
            continue
        z_range_global[angle_type] = (sub["depth_z"].min(), sub["depth_z"].max())
        th_range_global[angle_type] = (sub["theta"].min(), sub["theta"].max())

    # Pass 2: accumulate per-batch histograms
    hist_accum: dict = {}
    centers_cache: dict = {}

    for traj, df in traj_dfs.items():
        for f_start, f_end in BATCHES:
            rec = extract_records(df, frame_start=f_start, frame_end=f_end)
            if rec.empty:
                continue
            for angle_type in ("theta1", "theta2"):
                sub = rec[rec["angle_type"] == angle_type]
                if sub.empty:
                    continue
                z_r = z_range_global.get(angle_type)
                th_r = th_range_global.get(angle_type)
                hist, z_edges, th_edges = np.histogram2d(
                    sub["depth_z"], sub["theta"],
                    bins=[bins_z, bins_theta],
                    range=[z_r, th_r] if z_r and th_r else None,
                )
                if angle_type not in hist_accum:
                    hist_accum[angle_type] = []
                    centers_cache[angle_type] = (
                        (z_edges[:-1] + z_edges[1:]) / 2,
                        (th_edges[:-1] + th_edges[1:]) / 2,
                    )
                hist_accum[angle_type].append(hist)

    # Pass 3: mean + per-angle-type normalisation -> long DataFrame
    out_rows = []
    for angle_type, hists in hist_accum.items():
        mean_hist = np.mean(hists, axis=0)
        total = mean_hist.sum()
        if total > 0:
            mean_hist = mean_hist / total
        z_centers, th_centers = centers_cache[angle_type]
        for i, zc in enumerate(z_centers):
            for j, tc in enumerate(th_centers):
                out_rows.append({
                    "z_center": zc,
                    "theta_center": tc,
                    "count": mean_hist[i, j],
                    "angle_type": angle_type,
                })
    return pd.DataFrame(out_rows)


def main():
    args = parse_args()
    aa = args.aa

    traj_dfs = {}
    for traj in TRAJS:
        try:
            traj_dfs[traj] = load_traj(aa, traj)
        except FileNotFoundError as e:
            print(f"Warning: {e}")

    if not traj_dfs:
        raise SystemExit(f"No trajectory CSV found for analog '{aa}'.")

    freq = compute_batch_mean_histogram(traj_dfs, args.bins_z, args.bins_theta)

    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, f"freq_angle_{aa.lower()}.dat")
    freq.to_csv(out, sep="\t", index=False, float_format="%.10f")
    print(f"Frequencies written to '{out}' ({len(freq)} rows)")


if __name__ == "__main__":
    main()







