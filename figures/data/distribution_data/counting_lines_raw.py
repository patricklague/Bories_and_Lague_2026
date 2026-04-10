#!/usr/bin/env python3
"""
Scan raw_data/*/ and produce a summary .dat with:
  residu, traj, n_lines, n_frames
"""

import os
import glob
import numpy as np

RAW_DIR = "raw_data"
OUTPUT  = "raw_data_summary.dat"

rows = []

for resid_dir in sorted(glob.glob(os.path.join(RAW_DIR, "*"))):
    if not os.path.isdir(resid_dir):
        continue
    residu = os.path.basename(resid_dir)

    for dat_file in sorted(glob.glob(os.path.join(resid_dir, f"{residu}_contacts_*.dat"))):
        # Extract trajectory number from filename
        fname = os.path.basename(dat_file)
        traj = fname.replace(f"{residu}_contacts_", "").replace(".dat", "")

        data = np.loadtxt(dat_file, skiprows=1)
        n_lines = data.shape[0]
        n_frames = int(data[:, 0].max()) + 1

        rows.append((residu, traj, n_lines, n_frames))
        print(f"{residu} traj {traj}: {n_lines} lines, {n_frames} frames")

with open(OUTPUT, "w") as f:
    f.write("residu\ttraj\tn_lines\tn_frames\n")
    for residu, traj, n_lines, n_frames in rows:
        f.write(f"{residu}\t{traj}\t{n_lines}\t{n_frames}\n")

print(f"\nWrote {len(rows)} entries to {OUTPUT}")
