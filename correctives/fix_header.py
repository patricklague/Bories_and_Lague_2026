#!/usr/bin/env python3
"""
fix_header.py
=============
Normalise the headers of every data file in the repository so that they all
follow the canonical conventions defined in the check plan:

  * `frames` column is always called `#frame` (never `#section`)
  * trajectory-indexed columns use a hyphen: `<name>-1`, `<name>-2`, `<name>-3`
  * sub-block columns use the time-window naming
        bloc1 -> 401-600ns
        bloc2 -> 601-800ns
        bloc3 -> 801-1000ns
  * z-distribution / pmf files use the canonical header
        z  401-600ns  601-800ns  801-1000ns

The script *only* rewrites the first line (header) of each .dat file. The
numerical content is left untouched. A backup of every modified file is
written next to it with the suffix `.bak`.

NOTE: This script operates on the already-extracted `.dat` files inside this
repository. The raw trajectories themselves (PSF + compressed `.dcd.gz`
sections of 100 ns each) are hosted on Borealis under `POPC-<AA>/` and are
not touched here. See `correctives/list_of_script.md` for the trajectory
layout and the `DIR=...` snippet to paste into the analysis scripts.

Run from the repo root:
    python correctives/fix_header.py
"""

from __future__ import annotations

import shutil
from pathlib import Path

# ---------------------------------------------------------------------------
# Repository root (the script lives in <root>/correctives/)
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent


# ===========================================================================
# 1) MEMBRANE THICKNESS
# ---------------------------------------------------------------------------
# Path of the data to be modified:
#   figures/data/SUPP_membrane_parm/thickness/<aa>-thickness.dat
# ===========================================================================

def fix_thickness_headers() -> None:
    folder = ROOT / "figures/data/SUPP_membrane_parm/thickness"
    canonical = "#frame\tthickness-1\tthickness-2\tthickness-3\n"
    for f in sorted(folder.glob("*-thickness.dat")):
        _replace_first_line(f, canonical)

# Replaces:
#   "#frame  thickness-1  thickness-2  thickness-3"  (any whitespace) -> tab-separated canonical
#   any legacy "#section" -> "#frame"


# ===========================================================================
# 2) AREA PER LIPID
# ---------------------------------------------------------------------------
# Path of the data to be modified:
#   figures/data/SUPP_membrane_parm/area_per_lipid/<aa>-apl.dat
# ===========================================================================

def fix_apl_headers() -> None:
    folder = ROOT / "figures/data/SUPP_membrane_parm/area_per_lipid"
    canonical = (
        "#frame\tx-1\ty-1\tz-1\tapl-1"
        "\tx-2\ty-2\tz-2\tapl-2"
        "\tx-3\ty-3\tz-3\tapl-3\n"
    )
    for f in sorted(folder.glob("*-apl.dat")):
        _replace_first_line(f, canonical)

# Replaces:
#   "#section  x1 y1 z1 apl1  x2 y2 z2 apl2  x3 y3 z3 apl3" -> tab-separated canonical
#   "#section" -> "#frame"
#   "<col>1/2/3" -> "<col>-1/-2/-3" (consistent hyphen suffix per trajectory)


# ===========================================================================
# 3) DENSITY PROFILES (per analog, per component)
# ---------------------------------------------------------------------------
# Path of the data to be modified:
#   figures/data/SUPP_membrane_parm/densityProfiles/<aa>-<component>.dat
#   figures/data/densityProfile-popc/popc-<component>.dat
# ===========================================================================

def fix_densityProfiles_headers() -> None:
    folders = [
        ROOT / "figures/data/SUPP_membrane_parm/densityProfiles",
        ROOT / "figures/data/densityProfile-popc",
    ]
    bloc_map = {
        "bloc1": "401-600",
        "bloc2": "601-800",
        "bloc3": "801-1000",
    }
    for folder in folders:
        for f in sorted(folder.glob("*.dat")):
            with f.open("r", encoding="utf-8") as fh:
                header = fh.readline()
                rest = fh.read()
            new_header = header
            for old, new in bloc_map.items():
                new_header = new_header.replace(old, new)
            if new_header != header:
                _backup(f)
                with f.open("w", encoding="utf-8") as fh:
                    fh.write(new_header)
                    fh.write(rest)

# Replace : "bloc1 bloc2 bloc3"
# by      : "401-600ns 601-800ns 801-1000ns"
# Concretely, in the header line only:
#   "dens_traj<t>_bloc1" -> "dens_traj<t>_401-600ns"
#   "dens_traj<t>_bloc2" -> "dens_traj<t>_601-800ns"
#   "dens_traj<t>_bloc3" -> "dens_traj<t>_801-1000ns"
# (columns "z", "dens_mean", "se" are kept identical)


# ===========================================================================
# 4) LIPID ORDER PARAMETERS
# ---------------------------------------------------------------------------
# Path of the data to be modified:
#   figures/data/SUPP_membrane_parm/order_parameter/<aa>-chain{2,3}.dat
# ===========================================================================

def fix_order_parameter_headers() -> None:
    folder = ROOT / "figures/data/SUPP_membrane_parm/order_parameter"
    bloc_map = {
        "bloc1": "401-600",
        "bloc2": "601-800",
        "bloc3": "801-1000",
    }
    for f in sorted(folder.glob("*-chain*.dat")):
        with f.open("r", encoding="utf-8") as fh:
            header = fh.readline()
            rest = fh.read()
        new_header = header
        for old, new in bloc_map.items():
            new_header = new_header.replace(old, new)
        if new_header != header:
            _backup(f)
            with f.open("w", encoding="utf-8") as fh:
                fh.write(new_header)
                fh.write(rest)

# Replace : "bloc1 bloc2 bloc3"
# by      : "401-600ns 601-800ns 801-1000ns"
# Concretely, in the header line only:
#   "SCD_traj<t>_bloc1" -> "SCD_traj<t>_401-600ns"
#   "SCD_traj<t>_bloc2" -> "SCD_traj<t>_601-800ns"
#   "SCD_traj<t>_bloc3" -> "SCD_traj<t>_801-1000ns"
# (the "SCD_" prefix is the deuterium order parameter, not the analog name; kept)


# ===========================================================================
# 5) DISTRIBUTION DATA (total / monomer_4.5A / multimer_4.5A)
# ---------------------------------------------------------------------------
# Path of the data to be modified:
#   figures/data/distribution_data/total/<aa>/trajectory{1,2,3}.dat
#   figures/data/distribution_data/monomer_4.5A/<aa>/trajectory{1,2,3}.dat
#   figures/data/distribution_data/multimer_4.5A/<aa>/trajectory{1,2,3}.dat
#   figures/data/distribution_data/<mode>/<aa>/summary_<aa>.dat
# ===========================================================================

def fix_distribution_headers() -> None:
    canonical_traj = "z\t401-600ns\t601-800ns\t801-1000ns\n"
    base = ROOT / "figures/data/distribution_data"
    for mode in ("total", "monomer_4.5A", "multimer_4.5A"):
        for aa_dir in sorted((base / mode).glob("*/")):
            for f in sorted(aa_dir.glob("trajectory*.dat")):
                _replace_first_line(f, canonical_traj)
            # summary_<aa>.dat: rewrite per-traj/per-batch columns
            for f in sorted(aa_dir.glob("summary_*.dat")):
                _rename_traj_bloc_in_header(f)

# Replace : "bloc1 bloc2 bloc3"
# by      : "401-600ns 601-800ns 801-1000ns"
# Concretely, in the header line only:
#   trajectory<t>.dat header -> "z  401-600ns  601-800ns  801-1000ns"
#   summary_<aa>.dat:
#       "*_traj<t>_bloc1" -> "*_traj<t>_401-600ns"
#       "*_traj<t>_bloc2" -> "*_traj<t>_601-800ns"
#       "*_traj<t>_bloc3" -> "*_traj<t>_801-1000ns"


# ===========================================================================
# 6) PMF DATA (total / monomer_4.5A / multimer_4.5A)
# ---------------------------------------------------------------------------
# Path of the data to be modified:
#   figures/data/pmf_data/total/<aa>/trajectory{1,2,3}.dat
#   figures/data/pmf_data/monomer_4.5A/<aa>/trajectory{1,2,3}.dat
#   figures/data/pmf_data/multimer_4.5A/<aa>/trajectory{1,2,3}.dat
#   figures/data/pmf_data/<mode>/<aa>/pmf_<aa>.dat
# ===========================================================================

def fix_pmf_headers() -> None:
    canonical_traj = "z\t401-600ns\t601-800ns\t801-1000ns\n"
    base = ROOT / "figures/data/pmf_data"
    for mode in ("total", "monomer_4.5A", "multimer_4.5A"):
        for aa_dir in sorted((base / mode).glob("*/")):
            for f in sorted(aa_dir.glob("trajectory*.dat")):
                _replace_first_line(f, canonical_traj)
            for f in sorted(aa_dir.glob("pmf_*.dat")):
                _rename_traj_bloc_in_header(f)

# Replace : "bloc1 bloc2 bloc3"
# by      : "401-600ns 601-800ns 801-1000ns"
# Concretely, in the header line only:
#   trajectory<t>.dat header -> "z  401-600ns  601-800ns  801-1000ns"
#   pmf_<aa>.dat:
#       "*_traj<t>_bloc1" -> "*_traj<t>_401-600ns"
#       "*_traj<t>_bloc2" -> "*_traj<t>_601-800ns"
#       "*_traj<t>_bloc3" -> "*_traj<t>_801-1000ns"


# ===========================================================================
# 7) RAW CONTACT DATA
# ---------------------------------------------------------------------------
# Path of the data to be modified:
#   figures/data/distribution_data/raw_data/<aa>/<aa>_contacts_{1,2,3}.dat
# ===========================================================================

def fix_raw_contacts_headers() -> None:
    folder = ROOT / "figures/data/distribution_data/raw_data"
    canonical = "#frame\tresid\tx\ty\tz\t4.5A_cutoff\n"
    # Original column layout (0-indexed) produced by the contacts-vmd analysis:
    #   0: frame   1: resid   2: x   3: y   4: z
    #   5: 4.5A_cutoff   6: 6A_cutoff   7: 8A_cutoff   8: 10A_cutoff
    # We keep only columns [0, 1, 2, 3, 4, 5]
    # (frame, resid, x, y, z, 4.5A_cutoff) and drop the 6A / 8A / 10A
    # cutoff columns.
    keep_cols = (0, 1, 2, 3, 4, 5)
    for aa_dir in sorted(folder.glob("*/")):
        for f in sorted(aa_dir.glob("*_contacts_*.dat")):
            _backup(f)
            with f.open("r", encoding="utf-8") as fh:
                _ = fh.readline()  # discard old header
                lines = fh.readlines()
            with f.open("w", encoding="utf-8") as fh:
                fh.write(canonical)
                for ln in lines:
                    parts = ln.split()
                    if len(parts) < max(keep_cols) + 1:
                        # Not enough columns to safely subset -> keep line as-is.
                        fh.write(ln)
                        continue
                    fh.write("\t".join(parts[i] for i in keep_cols) + "\n")

# Replaces:
#   header line  -> canonical "#frame  resid  x  y  z  4.5A_cutoff"
#                   (6 columns, tab-sep)
#   data lines   -> only columns [frame, resid, x, y, z, 4.5A_cutoff] are
#                   kept; the 6A_cutoff / 8A_cutoff / 10A_cutoff columns
#                   are removed.


# ===========================================================================
# 8) AROMATICS ORIENTATION
# ---------------------------------------------------------------------------
# Path of the data to be modified:
#   figures/data/aromatics_orientation/freq_angle_<aa>.dat
# ===========================================================================

def fix_aromatics_headers() -> None:
    folder = ROOT / "figures/data/aromatics_orientation"
    canonical = "z_center\ttheta_center\tcount\tangle_type\n"
    for f in sorted(folder.glob("freq_angle_*.dat")):
        _replace_first_line(f, canonical)

# Replaces:
#   any first-line variant -> "z_center  theta_center  count  angle_type"



# ===========================================================================
# 9) AROMATICS ATOM COORDINATES — RENAME SECTION SUFFIX
# ---------------------------------------------------------------------------
# Path of the data to be renamed:
#   figures/data/aromatics_orientation/raw_data/total/<aa>/atom_coordinates/
#       traj{1,2,3}/<atom>_coor_{3,4,5}.dat
# Renames to:
#   figures/data/aromatics_orientation/raw_data/total/<aa>/atom_coordinates/
#       traj{1,2,3}/<atom>_coor_{401-600,601-800,801-1000}.dat
# ===========================================================================

def fix_atom_coordinates_filenames() -> None:
    base = ROOT / "figures/data/aromatics_orientation/raw_data/total"
    section_map = {
        "3": "401-600",
        "4": "601-800",
        "5": "801-1000",
    }
    for aa_dir in sorted(base.glob("*/")):
        coord_root = aa_dir / "atom_coordinates"
        if not coord_root.is_dir():
            continue
        for traj_dir in sorted(coord_root.glob("traj*/")):
            for f in sorted(traj_dir.glob("*_coor_*.dat")):
                stem = f.stem  # e.g. "CG_coor_3"
                # Only rename if the suffix after the last "_" is one of 3/4/5
                prefix, _, suffix = stem.rpartition("_")
                if suffix not in section_map:
                    continue
                new_name = f"{prefix}_{section_map[suffix]}{f.suffix}"
                new_path = f.with_name(new_name)
                if new_path.exists():
                    # Skip rather than clobber an already-renamed file
                    continue
                f.rename(new_path)


# ===========================================================================
# Helpers
# ===========================================================================

def _backup(path: Path) -> None:
    bak = path.with_suffix(path.suffix + ".bak")
    if not bak.exists():
        shutil.copy2(path, bak)


def _replace_first_line(path: Path, new_first_line: str) -> None:
    """Replace the first line of `path` with `new_first_line` (must end with \\n)."""
    with path.open("r", encoding="utf-8") as fh:
        _ = fh.readline()
        rest = fh.read()
    _backup(path)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(new_first_line)
        fh.write(rest)


def _rename_traj_bloc_in_header(path: Path) -> None:
    """Replace bloc1/2/3 by ns ranges in the header line only."""
    bloc_map = {
        "bloc1": "401-600ns",
        "bloc2": "601-800ns",
        "bloc3": "801-1000ns",
    }
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline()
        rest = fh.read()
    new_header = header
    for old, new in bloc_map.items():
        new_header = new_header.replace(old, new)
    if new_header != header:
        _backup(path)
        with path.open("w", encoding="utf-8") as fh:
            fh.write(new_header)
            fh.write(rest)


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    fix_thickness_headers()
    fix_apl_headers()
    fix_densityProfiles_headers()
    fix_order_parameter_headers()
    fix_distribution_headers()
    fix_pmf_headers()
    fix_raw_contacts_headers()
    fix_aromatics_headers()
    fix_atom_coordinates_filenames()
    print("Header normalisation complete. Backups written next to each modified file (*.bak).")


if __name__ == "__main__":
    main()