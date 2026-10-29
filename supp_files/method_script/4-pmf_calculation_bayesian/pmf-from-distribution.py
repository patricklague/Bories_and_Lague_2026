"""
PMF from the z-distributions with Bayesian (Dirichlet-multinomial) zero handling.
See PMF_bayesian_zero_handling.md for the equations and parameter choices.

Per analog:
1) Untranslated PMF of each 200 ns block (3 blocks x 3 trajectories).
   The stored probability histograms are converted back to exact counts using
   the number of molecule-frames per block (from the raw data). Instead of
   -kT ln(p + epsilon), each block uses the Dirichlet-multinomial
   posterior-mean proportion (Lovell, Chua & McGrath, NAR Genom. Bioinform.
   2, lqaa040, 2020):
       p_i = (n_i + alpha) / (N + K * alpha)
   with per-|z|-bin counts n_i, N = sum_i n_i, K bins and the Jeffreys
   pseudocount alpha = 0.5. Each block is treated on its own data only, so the
   block PMFs remain independent.
2) Block average and SE (std/sqrt(n)) of the PMFs without any translation.
3) Water reference a = inverse-variance weighted mean over [recal_lo, recal_hi].
4) Final PMF = average - a. SE (and sigma_a) from a block bootstrap repeating
   steps 2-4 on every replicate, x sqrt(n/(n-1)).
5) Outputs in <mode>/<analog>/ and <mode>/water_reference.dat.
"""
import os
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

# Distribution mode directory (e.g. "monomer_4.5A", "total", ...)
mode = "total"  # used to find distribution files and name output dir

# PMF constants
k_B = 0.008314  # kJ/(mol*K)
T = 303.15

# Dirichlet prior pseudocount: 0.5 = Jeffreys (default), 1.0 = Bayes-Laplace
bayes_alpha = 0.5
# Explicit sample size per block (molecule-frames); None = from the raw data
N_samples = None

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
data_root = os.path.join("..", "..", "..", "figures", "data", "distribution_data")
distribution_dir = os.path.join(data_root, mode)
raw_data_dir = os.path.join(data_root, "raw_data")
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


raw_summary = read_raw_summary(raw_summary_file)


def block_depth(analog, traj):
    """
    Number of molecule-frames in one block of one trajectory (the histogram
    total before normalisation). N_samples if set, otherwise taken from
    raw_data_summary.dat or, if absent, by counting the lines of the raw file.
    Valid for mode 'total', where every frame contains all solute molecules;
    filtered modes (monomer/multimer) have block-dependent totals that are
    not recoverable from the raw line count.
    """
    if N_samples is not None:
        return int(N_samples)
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

# ============================================================
# Dirichlet-multinomial posterior mean (Lovell et al. 2020)
# ============================================================

def counts_to_pmfs(counts, alpha=bayes_alpha):
    """
    counts: (..., K) folded |z| counts of one or more blocks.
    Posterior-mean proportions p_i = (n_i + alpha) / (N + K alpha) of a
    multinomial with a symmetric Dirichlet(alpha) prior, applied to each block
    independently. Returns the absolute PMFs -kT ln p.
    """
    k = counts.shape[-1]
    n_total = counts.sum(axis=-1, keepdims=True)
    p = (counts + alpha) / (n_total + k * alpha)
    return -k_B * T * np.log(p)
# ============================================================
# Block statistics
# ============================================================

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

    # --- Read all trajectory files, build common |z| grid ---
    traj_data = {}  # traj -> (z_array, density_matrix with n_batches cols)
    all_abs_z = set()

    for traj in range(1, n_trajectories + 1):
        fpath = os.path.join(dist_analog_dir, f"trajectory{traj}.dat")
        if not os.path.isfile(fpath):
            print(f"WARNING: {fpath} not found, skipping.")
            continue
        arr = np.loadtxt(fpath, skiprows=1)  # z, batch1, batch2, batch3
        traj_data[traj] = (arr[:, 0], arr[:, 1:])
        all_abs_z.update(np.abs(arr[:, 0]).tolist())

    if not traj_data:
        print(f"WARNING: No trajectory data for {analog}, skipping.")
        continue

    # Common |z| grid (positive only), extended to 50 A. Extended bins are
    # genuine zero-count observations (the solute never reached them).
    z_pos = np.array(sorted(all_abs_z))
    z_pos = z_pos[z_pos >= 0]
    if z_pos[-1] < 50:
        extra = np.arange(z_pos[-1] + z_bin, 50 + z_bin, z_bin)
        z_pos = np.concatenate([z_pos, extra])

    # --- Exact folded counts of every block ---
    block_counts, block_traj = [], []
    for traj, (z, densities) in traj_data.items():
        n_b = block_depth(analog, traj)
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

    # --- Step 1: untranslated block PMFs (Dirichlet posterior mean, per block) ---
    all_pmfs = counts_to_pmfs(block_counts)

    n_empty_win = int(np.sum(block_counts[:, mask_win] == 0))
    if n_empty_win == block_counts[:, mask_win].size:
        print(f"WARNING: {analog}: the water window is empty in every block; the water "
              f"reference is set by the prior floor, not by a measurement.")
    elif n_empty_win:
        print(f"WARNING: {analog}: {n_empty_win}/{block_counts[:, mask_win].size} block/bin "
              f"cells in the water window are empty; their PMF is the prior floor "
              f"and the water reference is sampling-limited.")

    # --- Step 2: block-averaged PMF and its SE, without any translation ---
    mean_raw, se_raw = block_mean_se(all_pmfs)
    if np.any(se_raw[mask_win] <= 1e-9):
        print(f"WARNING: {analog}: zero SE in the water window, unweighted mean used for a.")

    # --- Step 3: water reference a = inverse-variance weighted mean over the window ---
    a_ref = float(water_reference(mean_raw, se_raw, mask_win))
    # Analytical sigma_a, valid only if the window bins were independent (for information)
    se_win = se_raw[mask_win]
    sigma_a_indep = float(1.0 / np.sqrt(np.sum(1.0 / se_win ** 2))) if np.all(se_win > 1e-9) else float("nan")

    # --- Step 4: translation y'_i = y_i - a ---
    mean_pmf = mean_raw - a_ref

    # Uncertainty of y'_i by block bootstrap: resample the block PMFs with
    # replacement and repeat steps 2-4 on each replicate (each block PMF uses
    # only its own counts, so no refit is needed). This propagates sigma_a
    # with its full covariance with y_i (reduced error inside the window).
    rng = np.random.default_rng(bootstrap_seed)
    idx = rng.integers(0, n_blocks, size=(n_bootstrap, n_blocks))
    boot_mean, boot_se = block_mean_se(all_pmfs[idx])
    boot_a = water_reference(boot_mean, boot_se, mask_win)
    boot_shifted = boot_mean - boot_a[:, None]

    corr = np.sqrt(n_blocks / (n_blocks - 1)) if small_sample_correction else 1.0
    se_pmf = corr * boot_shifted.std(axis=0, ddof=1)
    sigma_a = float(corr * boot_a.std(ddof=1))

    print(f"{analog}: a = {a_ref:.4f} +/- {sigma_a:.4f} kJ/mol "
          f"(bootstrap; independent-bin estimate {sigma_a_indep:.4f})")
    reference_rows.append((analog, n_blocks, a_ref, sigma_a, sigma_a_indep))

    # --- Step 5: write outputs ---
    # Per-trajectory files: each block PMF shifted by the same common reference a,
    # so they share the water = 0 scale of the summary and keep their spread.
    header = f"{'z':>12s} {'400-600ns':>12s} {'600-800ns':>12s} {'800-1000ns':>12s}"
    for traj in traj_data:
        cols = [all_pmfs[i] - a_ref for i in range(n_blocks) if block_traj[i] == traj]
        output_file = os.path.join(out_analog_dir, f"trajectory{traj}.dat")
        np.savetxt(output_file, np.column_stack([z_pos] + cols),
                   fmt="%12.6f", header=header, comments="")
        print(f"Written {output_file}")

    header = f"{'z':>12s} {'mean':>12s} {'se':>12s}"

    # Final PMF: translated mean y'_i, bootstrap-propagated SE, and the number
    # of blocks with zero counts in each bin (n_empty = n_blocks: model-based value)
    n_empty = np.sum(block_counts == 0, axis=0)
    summary_file = os.path.join(out_analog_dir, f"pmf_{analog}.dat")
    np.savetxt(summary_file, np.column_stack([z_pos, mean_pmf, se_pmf, n_empty]),
               fmt=["%12.6f", "%12.6f", "%12.6f", "%12d"],
               header=header + f" {'n_empty':>12s}", comments="")
    print(f"Written {summary_file}")

    # Untranslated block-averaged PMF and its initial SE (step 2), for traceability
    raw_file = os.path.join(out_analog_dir, f"pmf_{analog}_untranslated.dat")
    np.savetxt(raw_file, np.column_stack([z_pos, mean_raw, se_raw]),
               fmt="%12.6f", header=header, comments="")
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