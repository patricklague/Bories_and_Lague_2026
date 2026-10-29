"""
PMF from the z-distributions by pooled counting: all blocks are summed before the logarithm.

Per analog:
1) Exact folded |z| counts of each 200 ns block (3 blocks x 3 trajectories). The stored
   probability histograms are converted back to integer counts with the number of
   molecule-frames per block (raw_data_summary.dat, otherwise the raw line count).
2) Pooled counts n(z) = sum_b n_b(z): every recorded event is used once, including the
   blocks with empty bins, and no pseudocount is added. A bin with n(z) = 0 has an
   undefined PMF (written as nan): it was never visited.
3) Water reference: the bulk is flat, so the 40-50 A bins are pooled into one reference
   probability, a = -kT ln(mean_{z in window} P(z)). It is defined as soon as one block
   visited the window.
4) Final PMF = -kT ln P(z) - a = -kT ln[n(z) / mean_window n].
5) SE by block bootstrap: the 9 count histograms are resampled with replacement and
   steps 2-4 are repeated; replicates with no water count (or no count in a bin) are
   undefined there and excluded, and the fraction of defined replicates is reported.
   SE is multiplied by sqrt(n/(n-1)).
Outputs in <mode>/<analog>/ and <mode>/water_reference.dat.
"""
import os
import warnings
import numpy as np

# ============================================================
# USER PARAMETERS
# ============================================================

analogs = [
    "sca", "scv", "scl", "sci", "scc", "scm", "scs", "sct", "scq",
    "scn", "scf", "scy", "scw", "scp", "glyd", "schd", "sche", "scdn",
    "scen", "sckn", "scrn", "schp", "sccm", "scym", "scd", "sce", "sck", "scr",
    "scrn-1", "scw-1", "scy-1", "scm-1", "scf-1", "scl-1", "scv-1", "sci-1"
]
n_trajectories = 3
n_batches = 3  # columns in each trajectory file (400-600, 600-800, 800-1000 ns)

# Distribution mode directory: "total" or "total_cyl" (every frame contains all solute
# molecules, so the block depth is known from the raw data).
mode = "total"

# PMF constants
k_B = 0.008314  # kJ/(mol*K)
T = 303.15

# Bulk-water reference window (angstroms): PMF is set to 0 here
recal_lo = 40.0
recal_hi = 50.0

# Z-bin width (for extending the grid)
z_bin = 1.0

# Block bootstrap
n_bootstrap = 10000
bootstrap_seed = 2025
small_sample_correction = True

# Paths
data_root = os.path.join("..", "..", "..", "figures", "data", "distribution_data")
raw_summary_file = os.path.join(data_root, "raw_data_summary.dat")

# ============================================================
# Counts and sampling depth
# ============================================================

def read_raw_summary(path):
    """residu, traj -> number of molecule-frames (rows) in the raw data."""
    table = {}
    if os.path.isfile(path):
        with open(path) as fh:
            next(fh)
            for line in fh:
                parts = line.split()
                if len(parts) >= 3:
                    table[(parts[0], int(parts[1]))] = int(parts[2])
    return table


def block_depth(analog, traj, raw_summary, raw_data_dir):
    """Number of molecule-frames in one 200 ns block (histogram total before normalisation)."""
    n_lines = raw_summary.get((analog, traj))
    if n_lines is None:
        raw = os.path.join(raw_data_dir, analog, f"{analog}_contacts_{traj}.dat")
        if not os.path.isfile(raw):
            raise FileNotFoundError(f"cannot determine sample size: {raw} not found")
        with open(raw, "rb") as fh:
            n_lines = sum(1 for _ in fh) - 1
    if n_lines % n_batches:
        raise ValueError(f"{analog} traj {traj}: {n_lines} rows not divisible into {n_batches} blocks")
    return n_lines // n_batches


def density_to_counts(density, n_total, label):
    """Recover the integer histogram counts from a normalised histogram."""
    raw = np.asarray(density, dtype=float) * n_total
    counts = np.rint(raw)
    resid = np.max(np.abs(raw - counts))
    if resid > 0.45 or abs(counts.sum() - n_total) > 0.01 * n_total:
        print(f"WARNING: {label}: counts not recovered exactly "
              f"(max rounding residual {resid:.2f}, sum {counts.sum():.0f} vs N {n_total}).")
    return counts


def fold_counts(z, counts, z_pos):
    """Sum the +z and -z counts of the same |z| bin onto the z_pos grid."""
    sym = np.zeros(len(z_pos))
    idx = np.array([np.argmin(np.abs(z_pos - abs(zi))) for zi in z])
    np.add.at(sym, idx, counts)
    return sym


def same_as_total(fpath, mode):
    """True if a non-'total' distribution file is numerically identical to the 'total' one."""
    ref = fpath.replace(os.sep + mode + os.sep, os.sep + "total" + os.sep)
    if ref == fpath or not os.path.isfile(ref):
        return False
    a, b = np.loadtxt(fpath, skiprows=1), np.loadtxt(ref, skiprows=1)
    return a.shape == b.shape and np.array_equal(a, b)

# ============================================================
# Pooled-count PMF
# ============================================================

def pooled_pmf(counts, mask_win):
    """
    counts: (..., n_blocks, K). Returns the untranslated PMF -kT ln P (P = pooled
    histogram normalised to 1), the water reference a, and the translated PMF.
    Undefined values (no count) are nan.
    """
    pooled = counts.sum(axis=-2)
    total = pooled.sum(axis=-1, keepdims=True)
    water = pooled[..., mask_win].mean(axis=-1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        untranslated = np.where(pooled > 0, -k_B * T * np.log(pooled / total), np.nan)
        a = np.where(water > 0, -k_B * T * np.log(water / total), np.nan)
    return untranslated, a[..., 0], untranslated - a


def bootstrap_se(values, n_blocks):
    """SE over bootstrap replicates (axis 0), ignoring undefined replicates."""
    corr = np.sqrt(n_blocks / (n_blocks - 1)) if small_sample_correction else 1.0
    defined = np.isfinite(values)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-nan columns stay nan
        se = corr * np.nanstd(values, axis=0, ddof=1)
    return np.where(defined.sum(axis=0) > 1, se, np.nan), defined.mean(axis=0)

# ============================================================
# Processing
# ============================================================

def run(mode):
    distribution_dir = os.path.join(data_root, mode)
    raw_data_dir = os.path.join(data_root, "raw_data_cyl" if mode.endswith("_cyl") else "raw_data")
    raw_summary = read_raw_summary(raw_summary_file)
    reference_rows = []

    for analog in analogs:
        dist_analog_dir = os.path.join(distribution_dir, analog)
        out_analog_dir = os.path.join(mode, analog)
        os.makedirs(out_analog_dir, exist_ok=True)

        # --- Read all trajectory files, build common |z| grid ---
        traj_data, all_abs_z = {}, set()
        for traj in range(1, n_trajectories + 1):
            fpath = os.path.join(dist_analog_dir, f"trajectory{traj}.dat")
            if not os.path.isfile(fpath):
                print(f"WARNING: {fpath} not found, skipping.")
                continue
            if mode != "total" and same_as_total(fpath, mode):
                print(f"WARNING: {fpath} is identical to the 'total' file; "
                      f"it was probably not regenerated from the {mode} raw data.")
            arr = np.loadtxt(fpath, skiprows=1)  # z, batch1, batch2, batch3
            traj_data[traj] = (arr[:, 0], arr[:, 1:])
            all_abs_z.update(np.abs(arr[:, 0]).tolist())
        if not traj_data:
            print(f"WARNING: No trajectory data for {analog}, skipping.")
            continue

        # Common |z| grid, extended to 50 A with genuine zero-count bins
        z_pos = np.array(sorted(all_abs_z))
        z_pos = z_pos[z_pos >= 0]
        if z_pos[-1] < 50:
            z_pos = np.concatenate([z_pos, np.arange(z_pos[-1] + z_bin, 50 + z_bin, z_bin)])

        # --- Step 1: exact folded counts of every block ---
        block_counts, block_traj = [], []
        for traj, (z, densities) in traj_data.items():
            n_b = block_depth(analog, traj, raw_summary, raw_data_dir)
            for b in range(n_batches):
                c = density_to_counts(densities[:, b], n_b, f"{analog} traj {traj} block {b + 1}")
                block_counts.append(fold_counts(z, c, z_pos))
                block_traj.append(traj)
        block_counts = np.array(block_counts)
        n_blocks = block_counts.shape[0]
        if n_blocks < 2:
            print(f"WARNING: {analog} has fewer than 2 blocks, cannot estimate SE, skipping.")
            continue
        mask_win = (z_pos >= recal_lo) & (z_pos <= recal_hi)
        if not mask_win.any():
            raise ValueError(f"{analog}: no grid point in the [{recal_lo}, {recal_hi}] A window.")

        # --- Steps 2-4: pooled counts, water reference, translation ---
        untranslated, a_ref, mean_pmf = pooled_pmf(block_counts, mask_win)
        a_ref = float(a_ref)
        water_counts = block_counts[:, mask_win].sum(axis=1)
        n_water_blocks = int((water_counts > 0).sum())
        if n_water_blocks == 0:
            print(f"WARNING: {analog}: no count in the water window in any block; "
                  f"the PMF cannot be referenced (written as nan).")
        elif n_water_blocks < n_blocks:
            print(f"WARNING: {analog}: water window visited in {n_water_blocks}/{n_blocks} blocks.")

        # --- Step 5: block bootstrap ---
        rng = np.random.default_rng(bootstrap_seed)
        idx = rng.integers(0, n_blocks, size=(n_bootstrap, n_blocks))
        boot_untranslated, boot_a, boot_pmf = pooled_pmf(block_counts[idx], mask_win)
        se_pmf, boot_defined = bootstrap_se(boot_pmf, n_blocks)
        se_untranslated, _ = bootstrap_se(boot_untranslated, n_blocks)
        sigma_a, a_defined = bootstrap_se(boot_a[:, None], n_blocks)
        sigma_a, a_defined = float(sigma_a[0]), float(a_defined[0])

        print(f"{analog}: a = {a_ref:.4f} +/- {sigma_a:.4f} kJ/mol; water visited in "
              f"{n_water_blocks}/{n_blocks} blocks; a defined in {100 * a_defined:.1f}% of replicates")
        reference_rows.append((analog, n_blocks, a_ref, sigma_a, n_water_blocks,
                               int(water_counts.sum()), a_defined))

        # --- Outputs ---
        # Per-trajectory files: each block PMF -kT ln p_b (no pseudocount, nan if empty),
        # shifted by the common reference a. Diagnostic only: the estimator is the pooled PMF.
        with np.errstate(divide="ignore", invalid="ignore"):
            p_blocks = block_counts / block_counts.sum(axis=1, keepdims=True)
            block_pmfs = np.where(block_counts > 0, -k_B * T * np.log(p_blocks), np.nan) - a_ref
        header = f"{'z':>12s} {'400-600ns':>12s} {'600-800ns':>12s} {'800-1000ns':>12s}"
        for traj in traj_data:
            cols = [block_pmfs[i] for i in range(n_blocks) if block_traj[i] == traj]
            output_file = os.path.join(out_analog_dir, f"trajectory{traj}.dat")
            np.savetxt(output_file, np.column_stack([z_pos] + cols), fmt="%12.6f", header=header, comments="")

        # Final PMF, bootstrap SE, number of empty blocks per bin and fraction of
        # bootstrap replicates in which the bin is defined.
        n_empty = np.sum(block_counts == 0, axis=0)
        summary_file = os.path.join(out_analog_dir, f"pmf_{analog}.dat")
        np.savetxt(summary_file, np.column_stack([z_pos, mean_pmf, se_pmf, n_empty, boot_defined]),
                   fmt=["%12.6f", "%12.6f", "%12.6f", "%12d", "%12.4f"],
                   header=f"{'z':>12s} {'mean':>12s} {'se':>12s} {'n_empty':>12s} {'boot_defined':>12s}",
                   comments="")

        # Untranslated pooled PMF -kT ln P and its bootstrap SE, for traceability
        raw_file = os.path.join(out_analog_dir, f"pmf_{analog}_untranslated.dat")
        np.savetxt(raw_file, np.column_stack([z_pos, untranslated, se_untranslated]),
                   fmt="%12.6f", header=f"{'z':>12s} {'mean':>12s} {'se':>12s}", comments="")
        print(f"Written {summary_file}, {raw_file} and trajectory files")

    if reference_rows:
        reference_file = os.path.join(mode, "water_reference.dat")
        with open(reference_file, "w") as fh:
            fh.write(f"{'analog':>12s} {'n_blocks':>9s} {'a':>12s} {'sigma_a':>12s} "
                     f"{'n_water_blocks':>15s} {'water_counts':>13s} {'a_defined':>10s}\n")
            for name, nb, a, sa, nw, wc, ad in reference_rows:
                fh.write(f"{name:>12s} {nb:>9d} {a:12.6f} {sa:12.6f} {nw:>15d} {wc:>13d} {ad:10.4f}\n")
        print(f"Written {reference_file}")


if __name__ == "__main__":
    run(mode)
    print("Done.")
