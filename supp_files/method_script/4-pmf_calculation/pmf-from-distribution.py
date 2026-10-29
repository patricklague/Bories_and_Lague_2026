import numpy as np
import os

# ============================================================
# USER PARAMETERS
# ============================================================

analogs = [
    "sca", "scv", "scl", "sci", "scc", "scm", "scs", "sct", "scq",
    "scn", "scf", "scy", "scw", "scp", "glyd", "schd", "sche", "scdn",
    "scen", "sckn", "scrn", "schp", "sccm", "scym", "scd", "sce", "sck", "scr",
    "scrn-1", "scw-1", "scl-1", "scv-1", "sci-1", "scm-1", "scf-1", "scy-1"
]
n_trajectories = 3
n_batches = 3  # columns in each trajectory file (400-600, 600-800, 800-1000 ns)

# Distribution mode directory (e.g. "monomer_4.5A", "total", ...)
mode = "total"  # used to find distribution files and name output dir

# PMF constants
k_B = 0.008314  # kJ/(mol*K)
T = 303.15
epsilon = 1e-5

# Bulk-water reference window (angstroms): PMF is set to 0 here
recal_lo = 40.0
recal_hi = 50.0

# Z-bin width (for extending the grid)
z_bin = 1.0

# Block bootstrap (uncertainty propagation of the water reference)
n_bootstrap = 10000
bootstrap_seed = 2025
# Rescale bootstrap SEs by sqrt(n/(n-1)) so that, for a plain mean, they equal
# the usual unbiased SE s/sqrt(n) (removes the small-sample bias of the bootstrap).
small_sample_correction = True

# Paths
distribution_dir = os.path.join("..", "..", "..", "figures", "data", "distribution_data", mode)
batch_labels = ["400-600ns", "600-800ns", "800-1000ns"]

# ============================================================
# Helper: density -> symmetrised PMF on |z| grid
# ============================================================

def density_to_pmf(z, density, z_pos):
    """
    Symmetrise density(z) around z=0 onto |z| grid and compute the absolute
    PMF = -kT ln(p + epsilon). No translation is applied here: the water
    reference is subtracted later, once, on the block-averaged PMF.
    Returns pmf array aligned with z_pos.
    """
    # Map density onto |z| grid: average neg and pos sides
    sym_density = np.zeros_like(z_pos, dtype=float)
    counts = np.zeros_like(z_pos, dtype=float)

    for zi, di in zip(z, density):
        az = abs(zi)
        idx = np.argmin(np.abs(z_pos - az))
        sym_density[idx] += di
        counts[idx] += 1.0

    # Average where we have contributions from both sides
    nonzero = counts > 0
    sym_density[nonzero] /= counts[nonzero]

    # Re-normalise so sum = 1 on the half-grid
    total = sym_density.sum()
    if total > 0:
        sym_density /= total

    # PMF = -kT ln(p + epsilon)
    return -k_B * T * np.log(sym_density + epsilon)


def block_mean_se(pmfs):
    """Mean and standard error over blocks (axis -2). pmfs: (..., n_blocks, n_z)."""
    n = pmfs.shape[-2]
    return pmfs.mean(axis=-2), pmfs.std(axis=-2, ddof=1) / np.sqrt(n)


def water_reference(mean, se, mask_win, se_floor=1e-9):
    """
    Inverse-variance weighted mean of the PMF over the water window.
    mean, se: (..., n_z). Returns a with shape (...).
    Falls back to an unweighted mean if any SE in the window is <= se_floor
    (identical values across blocks, e.g. a degenerate bootstrap replicate or
    a window empty in every block; such SEs are floating-point noise).
    """
    m = mean[..., mask_win]
    s = se[..., mask_win]
    degenerate = np.any(s <= se_floor, axis=-1, keepdims=True)
    w = np.where(degenerate, 1.0, 1.0 / np.where(s > se_floor, s, 1.0) ** 2)
    return np.sum(w * m, axis=-1) / np.sum(w, axis=-1)

# ============================================================
# Processing
# ============================================================

reference_rows = []  # (analog, n_blocks, a, sigma_a, sigma_a_indep)

for analog in analogs:
    dist_analog_dir = os.path.join(distribution_dir, analog)
    out_analog_dir = os.path.join(mode, analog)
    os.makedirs(out_analog_dir, exist_ok=True)

    # --- Pass 1: read all trajectory files, build common |z| grid ---
    traj_data = {}  # traj -> (z_array, density_matrix with n_batches cols)
    all_abs_z = set()

    for traj in range(1, n_trajectories + 1):
        fpath = os.path.join(dist_analog_dir, f"trajectory{traj}.dat")
        if not os.path.isfile(fpath):
            print(f"WARNING: {fpath} not found, skipping.")
            continue
        arr = np.loadtxt(fpath, skiprows=1)  # z, batch1, batch2, batch3
        z = arr[:, 0]
        densities = arr[:, 1:]  # shape (n_z, 3)
        traj_data[traj] = (z, densities)
        all_abs_z.update(np.abs(z).tolist())

    if not traj_data:
        print(f"WARNING: No trajectory data for {analog}, skipping.")
        continue

    # Common |z| grid (positive only), extended to 50 A
    z_pos = np.array(sorted(all_abs_z))
    z_pos = z_pos[z_pos >= 0]
    if z_pos[-1] < 50:
        extra = np.arange(z_pos[-1] + z_bin, 50 + z_bin, z_bin)
        z_pos = np.concatenate([z_pos, extra])
    n_z = len(z_pos)

    # --- Step 1: absolute (untranslated) PMF of every 200 ns block ---
    all_pmfs = []        # all blocks, shape (n_blocks, n_z)
    traj_pmfs_map = {}   # traj -> list of its block PMFs

    for traj in range(1, n_trajectories + 1):
        if traj not in traj_data:
            continue
        z, densities = traj_data[traj]
        traj_pmfs = [density_to_pmf(z, densities[:, b], z_pos) for b in range(n_batches)]
        traj_pmfs_map[traj] = traj_pmfs
        all_pmfs.extend(traj_pmfs)

    all_pmfs = np.array(all_pmfs)
    n_blocks = all_pmfs.shape[0]
    if n_blocks < 2:
        print(f"WARNING: {analog} has fewer than 2 blocks, cannot estimate SE, skipping.")
        continue

    mask_win = (z_pos >= recal_lo) & (z_pos <= recal_hi)
    if not mask_win.any():
        raise ValueError(f"{analog}: no grid point in the [{recal_lo}, {recal_hi}] A window.")

    # --- Step 2: block-averaged PMF and its SE, without any translation ---
    mean_raw, se_raw = block_mean_se(all_pmfs)

    # Empty bins in the water window give PMF = -kT ln(epsilon): the reference then
    # depends on epsilon rather than on sampled data.
    pmf_floor = -k_B * T * np.log(epsilon)
    n_floor = int(np.sum(np.isclose(all_pmfs[:, mask_win], pmf_floor)))
    if n_floor:
        print(f"WARNING: {analog}: {n_floor}/{all_pmfs[:, mask_win].size} block/bin values in the "
              f"water window are empty (epsilon floor); the water reference is poorly sampled.")
    if np.any(se_raw[mask_win] <= 1e-9):
        print(f"WARNING: {analog}: zero SE in the water window, unweighted mean used for a.")

    # --- Step 3: water reference a = inverse-variance weighted mean over the window ---
    a_ref = water_reference(mean_raw, se_raw, mask_win)
    # Analytical sigma_a, valid only if the window bins were independent (for information)
    se_win = se_raw[mask_win]
    sigma_a_indep = float(1.0 / np.sqrt(np.sum(1.0 / se_win ** 2))) if np.all(se_win > 1e-9) else float("nan")

    # --- Step 4: translation y'_i = y_i - a ---
    mean_pmf = mean_raw - a_ref

    # Uncertainty of y'_i by block bootstrap: resample the n_blocks blocks with
    # replacement and repeat steps 2-4 on each replicate. This propagates sigma_a
    # with the full covariance between y_i and a (reduced error inside the window,
    # where y_i and a are correlated; added error elsewhere).
    rng = np.random.default_rng(bootstrap_seed)
    idx = rng.integers(0, n_blocks, size=(n_bootstrap, n_blocks))
    boot_mean, boot_se = block_mean_se(all_pmfs[idx])     # (n_bootstrap, n_z)
    boot_a = water_reference(boot_mean, boot_se, mask_win)  # (n_bootstrap,)
    boot_shifted = boot_mean - boot_a[:, None]

    corr = np.sqrt(n_blocks / (n_blocks - 1)) if small_sample_correction else 1.0
    se_pmf = corr * boot_shifted.std(axis=0, ddof=1)
    sigma_a = corr * boot_a.std(ddof=1)

    print(f"{analog}: a = {a_ref:.4f} +/- {sigma_a:.4f} kJ/mol "
          f"(bootstrap; independent-bin estimate {sigma_a_indep:.4f})")
    reference_rows.append((analog, n_blocks, a_ref, sigma_a, sigma_a_indep))

    # --- Step 5: write outputs ---
    # Per-trajectory files: each block PMF shifted by the same common reference a,
    # so they share the water = 0 scale of the summary and keep their spread.
    header = f"{'z':>12s} {'400-600ns':>12s} {'600-800ns':>12s} {'800-1000ns':>12s}"
    for traj, traj_pmfs in traj_pmfs_map.items():
        output_file = os.path.join(out_analog_dir, f"trajectory{traj}.dat")
        np.savetxt(
            output_file,
            np.column_stack([z_pos] + [p - a_ref for p in traj_pmfs]),
            fmt="%12.6f",
            header=header,
            comments="",
        )
        print(f"Written {output_file}")

    header = f"{'z':>12s} {'mean':>12s} {'se':>12s}"

    # Final PMF: translated mean y'_i and bootstrap-propagated SE
    summary_file = os.path.join(out_analog_dir, f"pmf_{analog}.dat")
    np.savetxt(
        summary_file,
        np.column_stack([z_pos, mean_pmf, se_pmf]),
        fmt="%12.6f",
        header=header,
        comments="",
    )
    print(f"Written {summary_file}")

    # Untranslated block-averaged PMF and its initial SE (step 2), for traceability
    raw_file = os.path.join(out_analog_dir, f"pmf_{analog}_untranslated.dat")
    np.savetxt(
        raw_file,
        np.column_stack([z_pos, mean_raw, se_raw]),
        fmt="%12.6f",
        header=header,
        comments="",
    )
    print(f"Written {raw_file}")

# Water reference values for all analogs
if reference_rows:
    reference_file = os.path.join(mode, "water_reference.dat")
    with open(reference_file, "w") as fh:
        fh.write(f"{'analog':>12s} {'n_blocks':>9s} {'a':>12s} {'sigma_a':>12s} {'sigma_a_indep':>14s}\n")
        for name, nb, a, sa, sai in reference_rows:
            fh.write(f"{name:>12s} {nb:>9d} {a:12.6f} {sa:12.6f} {sai:14.6f}\n")
    print(f"Written {reference_file}")

print("Done.")