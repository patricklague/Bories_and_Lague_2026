"""
PMF from the z-distributions with the latent logarithm (lag) for empty bins.

Per analog:
1) Untranslated PMF of each 200 ns block (3 blocks x 3 trajectories).
   The stored probability histograms are converted back to exact counts using
   the number of molecule-frames per block (from the raw data). Instead of
   -kT ln(p + epsilon), the log-probability of every |z| bin is the latent
   logarithm of Biswas (2016, arXiv:1605.06064, "The latent logarithm"):
       z_bk ~ N(mu, sigma^2)                 (log latent rate)
       t_bk ~ Poisson(exp(z_bk) * o_b)        (observed count, depth o_b)
   fitted by iterated conditional modes (Newton MAP for z, sample mean and
   variance for mu, sigma^2). As recommended by Biswas for linked datasets,
   the 9 blocks share one latent prior (lag applied to their concatenation).
   Populated bins give lag ~ log(t/o); empty bins get a finite value that
   depends on the sampling depth and on the learned prior.
2) Block average and SE (std/sqrt(n)) of the PMFs without any translation.
3) Water reference a = inverse-variance weighted mean over [recal_lo, recal_hi].
4) Final PMF = average - a. SE (and sigma_a) from a block bootstrap repeating
   steps 1-4 (including the lag fit) on every replicate, x sqrt(n/(n-1)).
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

# Bulk-water reference window (angstroms): PMF is set to 0 here
recal_lo = 40.0
recal_hi = 50.0

# Z-bin width (for extending the grid)
z_bin = 1.0

# Latent logarithm (ICM) convergence
lag_tol = 1e-8
lag_max_iter = 5000

# Block bootstrap (uncertainty propagation of the water reference)
n_bootstrap = 10000
bootstrap_seed = 2025
# Rescale bootstrap SEs by sqrt(n/(n-1)) so that, for a plain mean, they equal
# the usual unbiased SE s/sqrt(n) (removes the small-sample bias of the bootstrap).
small_sample_correction = True
# Replicates processed together (memory/speed trade-off)
bootstrap_chunk = 1000

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
    total before normalisation). Taken from raw_data_summary.dat or, if
    absent, by counting the lines of the raw file. Valid for mode 'total',
    where every frame contains all solute molecules.
    """
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
# Latent logarithm (Biswas 2016)
# ============================================================

def _map_log_rate(t, o, mu, s2, z, tol=1e-10, max_newton=200, max_step=5.0):
    """
    MAP of z = log latent rate for t ~ Poisson(exp(z) o), z ~ N(mu, s2):
    maximise t z - o exp(z) - (z - mu)^2 / (2 s2) by Newton-Raphson
    (the objective is strictly concave, so the maximum is unique).
    """
    for _ in range(max_newton):
        e = o * np.exp(z)
        grad = t - e - (z - mu) / s2
        hess = -e - 1.0 / s2
        step = np.clip(grad / hess, -max_step, max_step)
        z = z - step
        if np.max(np.abs(step)) < tol:
            break
    return z


def latent_log(t, o, tol=lag_tol, max_iter=lag_max_iter):
    """
    Normalised latent logarithm nlag(t) = z for counts t with exposures o,
    with a common Normal prior N(mu, s2) on the last axis (leading axes are
    independent fits). Iterated conditional modes as in Biswas (2016),
    initialised with the sample mean and variance of log(t/o + 1).
    Exposures should be in units where the rates are O(1) for this init.
    Returns z (same shape as t), mu, s2 (shape t.shape[:-1] + (1,)).
    """
    t = np.asarray(t, dtype=float)
    o = np.broadcast_to(np.asarray(o, dtype=float), t.shape)
    y = np.log(t / o + 1.0)
    mu = y.mean(axis=-1, keepdims=True)
    s2 = np.maximum(y.var(axis=-1, keepdims=True), 1e-6)
    z = np.broadcast_to(mu, t.shape).copy()
    for _ in range(max_iter):
        z = _map_log_rate(t, o, mu, s2, z)
        mu_new = z.mean(axis=-1, keepdims=True)
        s2_new = z.var(axis=-1, keepdims=True)
        change = max(np.max(np.abs(mu_new - mu)), np.max(np.abs(s2_new - s2) / s2))
        mu, s2 = mu_new, s2_new
        if change < tol:
            break
    else:
        print(f"WARNING: latent logarithm not converged after {max_iter} ICM iterations "
              f"(last relative change {change:.2e}).")
    return z, mu, s2


def counts_to_pmfs(counts, depths):
    """
    counts: (..., n_blocks, K) folded |z| counts; depths: (..., n_blocks) N_b.
    Fits one lag prior over all blocks and bins of each leading index and
    returns absolute PMFs -kT ln p (p normalised to 1 per block), plus mu, s2.
    """
    k = counts.shape[-1]
    lead = counts.shape[:-2]
    # Exposure in units of the expected count per bin of a uniform histogram,
    # so latent rates are relative to uniform (O(1)), as in Biswas' examples.
    o = np.broadcast_to((depths / k)[..., None], counts.shape)
    z, mu, s2 = latent_log(counts.reshape(lead + (-1,)), o.reshape(lead + (-1,)))
    z = z.reshape(counts.shape)
    log_p = z - np.log(np.sum(np.exp(z), axis=-1, keepdims=True))
    return -k_B * T * log_p, mu, s2

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

reference_rows = []  # (analog, n_blocks, a, sigma_a, sigma_a_indep, mu, sigma)

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

    # --- Exact folded counts and sampling depth of every block ---
    block_counts, block_depths, block_traj = [], [], []
    for traj, (z, densities) in traj_data.items():
        n_b = block_depth(analog, traj)
        for b in range(n_batches):
            c = density_to_counts(densities[:, b], n_b, f"{analog} traj {traj} block {b + 1}")
            block_counts.append(fold_counts(z, c, z_pos))
            block_depths.append(c.sum())
            block_traj.append(traj)
    block_counts = np.array(block_counts)
    block_depths = np.array(block_depths, dtype=float)
    n_blocks = block_counts.shape[0]
    if n_blocks < 2:
        print(f"WARNING: {analog} has fewer than 2 blocks, cannot estimate SE, skipping.")
        continue

    mask_win = (z_pos >= recal_lo) & (z_pos <= recal_hi)
    if not mask_win.any():
        raise ValueError(f"{analog}: no grid point in the [{recal_lo}, {recal_hi}] A window.")

    # --- Step 1: untranslated block PMFs from the latent logarithm ---
    all_pmfs, lag_mu, lag_s2 = counts_to_pmfs(block_counts, block_depths)
    lag_mu, lag_sigma = float(lag_mu.squeeze()), float(np.sqrt(lag_s2.squeeze()))

    n_empty_win = int(np.sum(block_counts[:, mask_win] == 0))
    if n_empty_win == block_counts[:, mask_win].size:
        print(f"WARNING: {analog}: the water window is empty in every block; the water "
              f"reference is a pure latent-log extrapolation, not a measurement.")
    elif n_empty_win:
        print(f"WARNING: {analog}: {n_empty_win}/{block_counts[:, mask_win].size} block/bin "
              f"cells in the water window are empty; their PMF is a latent-log estimate "
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

    # Uncertainty of y'_i by block bootstrap: resample the blocks (raw counts)
    # with replacement and repeat steps 1-4 on each replicate, including the
    # lag fit whose shared prior couples the blocks. This propagates sigma_a
    # with its full covariance with y_i (reduced error inside the window).
    rng = np.random.default_rng(bootstrap_seed)
    idx = rng.integers(0, n_blocks, size=(n_bootstrap, n_blocks))
    boot_shifted = np.empty((n_bootstrap, len(z_pos)))
    boot_a = np.empty(n_bootstrap)
    for start in range(0, n_bootstrap, bootstrap_chunk):
        sl = slice(start, start + bootstrap_chunk)
        pmfs_b, _, _ = counts_to_pmfs(block_counts[idx[sl]], block_depths[idx[sl]])
        bm, bs = block_mean_se(pmfs_b)
        ba = water_reference(bm, bs, mask_win)
        boot_a[sl] = ba
        boot_shifted[sl] = bm - ba[:, None]

    corr = np.sqrt(n_blocks / (n_blocks - 1)) if small_sample_correction else 1.0
    se_pmf = corr * boot_shifted.std(axis=0, ddof=1)
    sigma_a = float(corr * boot_a.std(ddof=1))

    print(f"{analog}: a = {a_ref:.4f} +/- {sigma_a:.4f} kJ/mol "
          f"(bootstrap; independent-bin estimate {sigma_a_indep:.4f}); "
          f"lag prior mu = {lag_mu:.3f}, sigma = {lag_sigma:.3f}")
    reference_rows.append((analog, n_blocks, a_ref, sigma_a, sigma_a_indep, lag_mu, lag_sigma))

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

# Water reference values and lag prior for all analogs
if reference_rows:
    reference_file = os.path.join(mode, "water_reference.dat")
    with open(reference_file, "w") as fh:
        fh.write(f"{'analog':>12s} {'n_blocks':>9s} {'a':>12s} {'sigma_a':>12s} "
                 f"{'sigma_a_indep':>14s} {'lag_mu':>10s} {'lag_sigma':>10s}\n")
        for name, nb, a, sa, sai, mu, sg in reference_rows:
            fh.write(f"{name:>12s} {nb:>9d} {a:12.6f} {sa:12.6f} {sai:14.6f} {mu:10.4f} {sg:10.4f}\n")
    print(f"Written {reference_file}")

print("Done.")