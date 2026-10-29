"""
PMF from the z-distributions with a Poisson-Gamma model for empty bins.
Cylinder version: reads distribution_data/total_cyl, writes total_cyl/.

Per analog:
1) Untranslated PMF of each 200 ns block (3 blocks x 3 trajectories).
   The stored probability histograms are converted back to exact counts using
   the number of molecule-frames per block (from the raw data). Instead of
   -kT ln(p + epsilon), each block's probability in every |z| bin is the
   latent abundance of the Poisson-Gamma (PG) model of Jiang, Zhan & Wang
   (2023, JASA 118:792, "A flexible zero-inflated Poisson-Gamma model"):
       W_bk | U_bk ~ Poisson(pi_k M_b U_bk),   U_bk ~ Gamma(1/theta_k, theta_k)
   with |z| bin k <-> taxon, block b <-> repeated measurement, trajectory <->
   subject, M_b = block depth (molecule-frames), pi_k = long-term probability
   of the bin and theta_k = bin-specific dispersion (var(U) = theta_k) that
   absorbs the block-to-block over-dispersion of correlated MD frames.
   The zero-inflation probability of ZIPG is set to 0: in an equilibrium
   Boltzmann distribution every bin has a non-zero probability, so empty bins
   arise only from finite sampling (no structural zeros); ZIPG then reduces
   to its PG part. pi_k and theta_k are maximum-likelihood estimates from the
   blocks; the block probability is the posterior mean of pi_k U_bk:
       p_bk = pi_k (W_bk + 1/theta_k) / (pi_k M_b + 1/theta_k)
   so an empty block is shrunk toward the bin's abundance in the other blocks.
   Bins empty in every block carry no information on pi_k; they receive the
   Jeffreys estimate pi_k = 0.5 / sum_b M_b (half a count over all blocks).
2) Block average of the PMFs without any translation. Because the PG fit
   pools the blocks and shrinks each toward the bin's pooled abundance, the
   block PMFs are not independent: their SE is the bootstrap SE of this
   untranslated average (not std/sqrt(n) across blocks).
3) Water reference a = inverse-variance weighted mean over [recal_lo, recal_hi],
   with weights 1/SE^2 from step 2, kept fixed in the bootstrap.
4) Final PMF = average - a. SE (and sigma_a) from a block bootstrap repeating
   steps 1-4 (including the PG fit) on every replicate, x sqrt(n/(n-1)).
5) Outputs in <mode>/<analog>/ and <mode>/water_reference.dat. The
   per-trajectory files contain the shrunk (posterior) block PMFs.
"""
import os
import numpy as np
from scipy.special import gammaln

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
mode = "total_cyl"  # used to find distribution files and name output dir

# PMF constants
k_B = 0.008314  # kJ/(mol*K)
T = 303.15

# Bulk-water reference window (angstroms): PMF is set to 0 here
recal_lo = 40.0
recal_hi = 50.0

# Z-bin width (for extending the grid)
z_bin = 1.0

# Poisson-Gamma dispersion fit (theta_k, maximum likelihood on log theta)
theta_min = 1e-8         # ~ Poisson limit (no over-dispersion)
theta_max = 1e4          # far above any fitted value (largest fitted ~ 46)
theta_log_tol = 1e-4     # golden-section tolerance on ln(theta)
pi_zero_pseudocount = 0.5  # Jeffreys estimate for bins empty in every block

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
raw_data_dir = os.path.join(data_root, "raw_data_cyl")
frames_per_block = 20000  # 200 ns blocks (10 ps per frame)
# ============================================================
# Counts and sampling depth
# ============================================================

def block_depth(analog, traj):
    """
    Number of molecule-frames in one block of one trajectory (the histogram
    total before normalisation) = molecules per frame x frames per block.
    Molecules per frame are read from the first frame of the cylinder raw
    file, so the value stays correct even if the raw file is incomplete; a
    warning is printed when its line count does not match. Valid for mode
    'total_cyl', where every frame contains all solute molecules.
    """
    raw = os.path.join(raw_data_dir, analog, f"{analog}_contacts_{traj}.dat")
    if not os.path.isfile(raw):
        raise FileNotFoundError(f"cannot determine sample size: {raw} not found")
    n_mol, n_lines, first = 0, 0, None
    with open(raw) as fh:
        next(fh)
        for line in fh:
            n_lines += 1
            frame = line.split(None, 1)[0]
            if first is None:
                first = frame
            if frame == first:
                n_mol += 1
    expected = n_mol * frames_per_block * n_batches
    if n_lines != expected:
        print(f"WARNING: {raw} has {n_lines} rows instead of {expected} "
              f"({n_mol} molecules x {frames_per_block * n_batches} frames); incomplete raw file.")
    return n_mol * frames_per_block

def same_as_total(fpath):
    """True if this total_cyl distribution is numerically identical to the non-cylinder one."""
    ref = fpath.replace(os.sep + mode + os.sep, os.sep + "total" + os.sep)
    if ref == fpath or not os.path.isfile(ref):
        return False
    a, b = np.loadtxt(fpath, skiprows=1), np.loadtxt(ref, skiprows=1)
    return a.shape == b.shape and np.array_equal(a, b)


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
# Poisson-Gamma model (Jiang, Zhan & Wang 2023, without zero inflation)
# ============================================================

def _nb_loglik(W, mu, log_theta):
    """
    Marginal log-likelihood of the PG model (negative binomial with mean mu and
    variance mu (1 + theta mu)), summed over blocks (axis -2), up to the
    theta-independent term -ln W!. W, mu: (..., n_blocks, K); log_theta: (..., 1, K).
    """
    theta = np.exp(log_theta)
    r = 1.0 / theta
    x = mu * theta
    ll = gammaln(W + r) - gammaln(r) - r * np.log1p(x) + W * (np.log(x) - np.log1p(x))
    return ll.sum(axis=-2, keepdims=True)


def fit_dispersion(W, mu):
    """
    Maximum-likelihood ln(theta_k) for every bin by golden-section search on
    [ln theta_min, ln theta_max] (vectorised over all leading axes and bins).
    Under-dispersed bins converge to theta_min (Poisson limit).
    """
    shape = W.shape[:-2] + (1, W.shape[-1])
    a = np.full(shape, np.log(theta_min))
    b = np.full(shape, np.log(theta_max))
    g = (np.sqrt(5.0) - 1.0) / 2.0
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = _nb_loglik(W, mu, c), _nb_loglik(W, mu, d)
    while np.max(b - a) > theta_log_tol:
        left = fc > fd                      # maximum lies in [a, d]
        b = np.where(left, d, b)
        a = np.where(left, a, c)
        c_new = np.where(left, b - g * (b - a), d)
        d_new = np.where(left, c, a + g * (b - a))
        f_new = _nb_loglik(W, mu, np.where(left, c_new, d_new))
        fc, fd = np.where(left, f_new, fd), np.where(left, fc, f_new)
        c, d = c_new, d_new
    return 0.5 * (a + b)


def counts_to_pmfs(counts, depths):
    """
    counts: (..., n_blocks, K) folded |z| counts; depths: (..., n_blocks) M_b.
    Fits the PG model bin by bin across the blocks of each leading index and
    returns absolute PMFs -kT ln p (p normalised to 1 per block), pi and theta.
    """
    M = depths[..., None]                                  # (..., n_blocks, 1)
    total = counts.sum(axis=-2, keepdims=True)             # (..., 1, K)
    # Long-term abundance: pooled ratio estimator (the ML estimate when all
    # blocks have the same depth, as here); Jeffreys value for all-empty bins.
    pi = np.where(total > 0, total, pi_zero_pseudocount) / M.sum(axis=-2, keepdims=True)
    mu = pi * M
    log_theta = fit_dispersion(counts, mu)
    log_theta = np.where(total > 0, log_theta, np.log(theta_min))
    r = np.exp(-log_theta)
    # Posterior mean of the block's latent abundance pi_k U_bk (Gamma-Poisson conjugacy)
    p = pi * (counts + r) / (mu + r)
    p = p / p.sum(axis=-1, keepdims=True)
    return -k_B * T * np.log(p), pi, np.exp(log_theta)
# ============================================================
# Block statistics
# ============================================================

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

reference_rows = []  # (analog, n_blocks, a, sigma_a, sigma_a_indep, theta_win)

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
        if same_as_total(fpath):
            print(f"WARNING: {fpath} is identical to the non-cylinder 'total' file; "
                  f"it was probably not regenerated from the cylinder raw data.")
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

    # --- Step 1: untranslated block PMFs from the Poisson-Gamma model ---
    all_pmfs, pg_pi, pg_theta = counts_to_pmfs(block_counts, block_depths)
    theta_win = float(np.median(pg_theta[0, mask_win]))

    n_empty_win = int(np.sum(block_counts[:, mask_win] == 0))
    if n_empty_win == block_counts[:, mask_win].size:
        print(f"WARNING: {analog}: the water window is empty in every block; the water "
              f"reference is a model extrapolation (Jeffreys value), not a measurement.")
    elif n_empty_win:
        print(f"WARNING: {analog}: {n_empty_win}/{block_counts[:, mask_win].size} block/bin "
              f"cells in the water window are empty; their PMF is a Poisson-Gamma "
              f"posterior estimate and the water reference is sampling-limited.")

    # Block bootstrap: resample the blocks (raw counts) with replacement and
    # repeat the PG fit and the block average on each replicate. The PG fit
    # pools information across blocks (pi_k, theta_k) and shrinks each block
    # toward the bin's pooled abundance, so the block PMFs are not independent
    # and std/sqrt(n) across them would underestimate the SE; the bootstrap
    # SE of the untranslated mean is used instead.
    rng = np.random.default_rng(bootstrap_seed)
    idx = rng.integers(0, n_blocks, size=(n_bootstrap, n_blocks))
    boot_mean = np.empty((n_bootstrap, len(z_pos)))
    for start in range(0, n_bootstrap, bootstrap_chunk):
        sl = slice(start, start + bootstrap_chunk)
        pmfs_b, _, _ = counts_to_pmfs(block_counts[idx[sl]], block_depths[idx[sl]])
        boot_mean[sl] = pmfs_b.mean(axis=-2)
    corr = np.sqrt(n_blocks / (n_blocks - 1)) if small_sample_correction else 1.0

    # --- Step 2: block-averaged PMF (no translation) and its bootstrap SE ---
    mean_raw = all_pmfs.mean(axis=0)
    se_raw = corr * boot_mean.std(axis=0, ddof=1)
    if np.any(se_raw[mask_win] <= 1e-9):
        print(f"WARNING: {analog}: zero SE in the water window, unweighted mean used for a.")

    # --- Step 3: water reference a = inverse-variance weighted mean over the window ---
    # The weights 1/se_raw^2 are fixed and reused in every bootstrap replicate
    # (re-estimating them per replicate would require a nested bootstrap).
    a_ref = float(water_reference(mean_raw, se_raw, mask_win))
    # Analytical sigma_a, valid only if the window bins were independent (for information)
    se_win = se_raw[mask_win]
    sigma_a_indep = float(1.0 / np.sqrt(np.sum(1.0 / se_win ** 2))) if np.all(se_win > 1e-9) else float("nan")

    # --- Step 4: translation y'_i = y_i - a, SE from the same bootstrap ---
    # Subtracting each replicate's own reference propagates sigma_a with its
    # full covariance with y_i (reduced error inside the window).
    mean_pmf = mean_raw - a_ref
    boot_a = water_reference(boot_mean, np.broadcast_to(se_raw, boot_mean.shape), mask_win)
    se_pmf = corr * (boot_mean - boot_a[:, None]).std(axis=0, ddof=1)
    sigma_a = float(corr * boot_a.std(ddof=1))

    print(f"{analog}: a = {a_ref:.4f} +/- {sigma_a:.4f} kJ/mol "
          f"(bootstrap; independent-bin estimate {sigma_a_indep:.4f}); "
          f"median theta in water window = {theta_win:.3g}")
    reference_rows.append((analog, n_blocks, a_ref, sigma_a, sigma_a_indep, theta_win))

    # --- Step 5: write outputs ---
    # Per-trajectory files: each block PMF (PG posterior, shrunk toward the
    # bin's pooled abundance) shifted by the same common reference a.
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

# Water reference values and PG dispersion for all analogs
if reference_rows:
    reference_file = os.path.join(mode, "water_reference.dat")
    with open(reference_file, "w") as fh:
        fh.write(f"{'analog':>12s} {'n_blocks':>9s} {'a':>12s} {'sigma_a':>12s} "
                 f"{'sigma_a_indep':>14s} {'theta_win':>12s}\n")
        for name, nb, a, sa, sai, th in reference_rows:
            fh.write(f"{name:>12s} {nb:>9d} {a:12.6f} {sa:12.6f} {sai:14.6f} {th:12.4g}\n")
    print(f"Written {reference_file}")

print("Done.")