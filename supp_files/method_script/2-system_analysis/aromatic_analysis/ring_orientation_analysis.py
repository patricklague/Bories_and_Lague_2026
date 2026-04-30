#!/usr/bin/env python3
"""
Script to analyze the orientation of aromatic rings from coordinate files.

For each trajectory (1, 2, 3) of the residue `aa`:
  - reads the coordinate files per atom in
        figures/data/aromatics_orientation/raw_data/total/{aa}/atom_coordinates/traj{traj}/
    (one file per atom and per 200 ns section:
     *_coor_401-600.dat, *_coor_601-800.dat, *_coor_801-1000.dat)
  - merges the sections by adjusting the frame indices,
  - calculates for each frame / each residue:
      * the angle between the normal to the aromatic ring plane and the Z axis
      * the angle between the vector (atom_pair) and the Z axis
      * the Z depth of the ring centroid
  - writes a file
        figures/data/aromatics_orientation/raw_data/total/{aa}/vector_orientations/orientation_{aa}_traj{traj}.csv

One execution thus produces 3 CSVs (one per trajectory).
"""
import os
import glob
import re
import sys
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Section labels (200 ns each). The script treats them as ordered "parts" and
# rescales the frame index across them so that the resulting frame numbering
# is continuous over the whole 401-1000 ns window.
# ---------------------------------------------------------------------------
SECTIONS = ["401-600", "601-800", "801-1000"]

# Repository layout (relative to this script)
DATA_ROOT = os.path.join(
    "..", "..", "..", "..",  # method_script/2-system_analysis/aromatic_analysis -> repo root
    "figures", "data", "aromatics_orientation", "raw_data", "total",
)


def parse_coor_files(aa: str, traj: int) -> dict:
    """Read all *_coor_<section>.dat files for one (aa, traj) and return
    {frame: {resid: {atom: np.array([x, y, z])}}}.

    Frames from successive sections are offset so the resulting numbering is
    continuous across the three 200 ns blocks.
    """
    coords: dict = {}
    max_frame = -1
    base = os.path.join(DATA_ROOT, aa, "atom_coordinates", f"traj{traj}")

    for section in SECTIONS:
        files = sorted(glob.glob(os.path.join(base, f"*_coor_{section}.dat")))
        if not files:
            continue
        offset = max_frame + 1 if max_frame >= 0 else 0
        section_max = max_frame
        for fp in files:
            atom = os.path.basename(fp).split("_")[0]
            with open(fp) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    frame, _, idx, resid, x, y, z = line.split()
                    frame = int(frame) + offset
                    resid = int(resid)
                    coords.setdefault(frame, {}).setdefault(resid, {})[atom] = (
                        np.array([float(x), float(y), float(z)])
                    )
                    section_max = max(section_max, frame)
        max_frame = section_max

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
    cos_t = np.dot(v, axis) / (norm_v * norm_axis)
    cos_t = np.clip(cos_t, -1.0, 1.0)
    return np.degrees(np.arccos(cos_t))


def analyze_coords(coords, atom_pair, ring_atoms):
    rows = []
    a1, a2 = atom_pair
    for frame in sorted(coords):
        for resid, atoms in coords[frame].items():
            pts = [atoms[a] for a in ring_atoms if a in atoms]
            if len(pts) < 3:
                continue
            centroid, normal = ring_normal(pts)
            depth_z = centroid[2]
            axis = np.array([0, 0, 1]) if depth_z >= 0 else np.array([0, 0, -1])

            if a1 in atoms and a2 in atoms:
                u = unit_vector(atoms[a1], atoms[a2])
                theta2 = angle_between(u, axis)
            else:
                theta2 = np.nan
            theta1 = angle_between(normal, axis)

            rows.append({
                "frame": frame,
                "index": resid,
                f"{a1}-{a2}_angle_z": theta2,
                "normal_angle_z": theta1,
                "depth_z": depth_z,
            })
    return pd.DataFrame(rows)


def main(aa: str | None = None):
    if aa in ("SCY", "SCF"):
        atom_pair = ("CG", "CZ")
        ring_atoms = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    elif aa == "SCW":
        atom_pair = ("CZ3", "CE2")
        ring_atoms = ["CD2", "CZ2", "CZ3", "CD1", "CE2", "CH2", "CE3"]
    else:
        print(f"Erreur: acide aminé '{aa}' non reconnu. Valeurs acceptées: SCY, SCF, SCW")
        return

    out_dir = os.path.join(DATA_ROOT, aa, "vector_orientations")
    os.makedirs(out_dir, exist_ok=True)

    for traj in (1, 2, 3):
        coords = parse_coor_files(aa, traj)
        if not coords:
            print(f"No data found for {aa} traj{traj}, skipping this trajectory.")
            continue
        df_traj = analyze_coords(coords, atom_pair, ring_atoms)
        if df_traj.empty:
            print(f"No data extracted for {aa} traj{traj}, skipping this trajectory.")
            continue
        df_traj = df_traj.sort_values(["frame", "index"]).reset_index(drop=True)

        out_path = os.path.join(out_dir, f"orientation_{aa.lower()}_traj{traj}.csv")
        df_traj.to_csv(out_path, index=False)
        print(f"Output file: {out_path}")


if __name__ == "__main__":
    aa = sys.argv[1] if len(sys.argv) > 1 else None
    main(aa)






