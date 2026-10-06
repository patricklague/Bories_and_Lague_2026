import numpy as np
import os

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

# ------------------------------------------------------------------
# Bayesian handling of empty (zero-count) bins
# ------------------------------------------------------------------
# Instead of adding an ad-hoc epsilon inside the logarithm, zero counts are
# treated with a Bayesian (Dirichlet-multinomial) posterior-mean proportion
# (Lovell, Chua & McGrath, "Counts: an outstanding challenge for log-ratio
# analysis of compositional data in the molecular biosciences",
# NAR Genomics and Bioinformatics 2(2), lqaa040, 2020):
#
#       p_i = (n_i + alpha) / (N + K * alpha)
#
# with per-|z|-bin counts n_i, total N = sum_i n_i, K bins and a Dirichlet
# prior pseudocount alpha. Empty bins then receive a finite, sample-size-aware
# floor alpha / (N + K*alpha) instead of the constant -kT ln(epsilon).
bayes_alpha = 0.5   # Dirichlet prior pseudocount: 1.0 = Bayes-Laplace, 0.5 = Jeffreys (default)
N_samples = None    # explicit per-histogram sample size (counts); None = infer from the histogram quantisation

# PMF recalibration window (angstroms)
recal_lo = 40.0
recal_hi = 50.0

# Z-bin width (for extending the grid)
z_bin = 1.0

# Paths
#distribution_dir = os.path.join("..", "distribution_data", mode)
distribution_dir = os.path.join("..", "..", "..", "figures", "data", "distribution_data", mode)
batch_labels = ["400-600ns", "600-800ns", "800-1000ns"]

# ============================================================
# Helper: density -> symmetrised PMF on |z| grid
# ============================================================

def density_to_pmf(z, density, z_pos, alpha=bayes_alpha, N=N_samples):
    """
    Symmetrise density(z) around z=0 onto the |z| grid and convert to an absolute
    PMF (-kT ln p). No reference offset is applied here: every batch of a given
    sidechain is later anchored at a single common reference bin (the minimum of
    the mean PMF, i.e. the most-populated depth) and the bulk-water zero is set
    once as a rigid shift of the averaged curve (see the processing loop).

    Empty (zero-count) bins are handled with a Bayesian (Dirichlet-multinomial)
    posterior mean rather than an ad-hoc epsilon inside the logarithm
    (Lovell, Chua & McGrath, NAR Genomics Bioinformatics 2, lqaa040, 2020):

        p_i = (n_i + alpha) / (N + K * alpha)

    where n_i are the per-|z|-bin counts, N = sum_i n_i, K the number of bins
    and alpha the Dirichlet prior pseudocount. Empty bins thus receive a finite
    floor alpha / (N + K*alpha) that scales with the actual sample size.
    Returns pmf array aligned with z_pos.
    """
    density = np.asarray(density, dtype=float)
    positive = density > 0

    # Reconstruct integer per-bin counts from the stored probability histogram.
    # The smallest positive density corresponds to a single observation (~1/N)
    # and sets the quantisation used to recover counts; N may also be supplied.
    if N is not None:
        counts_signed = np.rint(density * float(N))
    elif positive.any():
        quantum = density[positive].min()
        counts_signed = np.rint(density / quantum)
    else:
        counts_signed = np.zeros_like(density)

    # Fold the +z and -z observations of the same |z| by SUMMING their counts.
    sym_counts = np.zeros_like(z_pos, dtype=float)
    for zi, ci in zip(z, counts_signed):
        idx = np.argmin(np.abs(z_pos - abs(zi)))
        sym_counts[idx] += ci

    n_total = sym_counts.sum()
    n_bins = len(z_pos)

    # Dirichlet-multinomial posterior-mean proportions (strictly positive -> finite PMF).
    if n_total > 0:
        p = (sym_counts + alpha) / (n_total + n_bins * alpha)
    else:
        p = np.full(n_bins, 1.0 / n_bins)

    # PMF = -kT ln(p)  (absolute; common reference applied in the processing loop)
    pmf_raw = -k_B * T * np.log(p)

    return pmf_raw

# ============================================================
# Processing
# ============================================================

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

    # --- Pass 2: compute the absolute PMF of every batch ---
    all_pmfs = []           # all 9 batch PMFs (absolute, -kT ln p)
    traj_pmfs_map = {}      # traj -> list of batch PMFs (absolute)

    for traj in range(1, n_trajectories + 1):
        if traj not in traj_data:
            continue
        z, densities = traj_data[traj]
        traj_pmfs = []
        for b in range(n_batches):
            pmf = density_to_pmf(z, densities[:, b], z_pos)
            traj_pmfs.append(pmf)
            all_pmfs.append(pmf)
        traj_pmfs_map[traj] = traj_pmfs

    all_pmfs = np.array(all_pmfs)  # shape (up to 9, n_z)

    # --- Common reference = z of minimum energy (= most-populated depth) ---
    # The reference bin is the argmin of the across-batch mean PMF and is applied
    # identically to all 9 replicas, so every replica passes through 0 at ref_z
    # and the error bars vanish there and grow outward.
    mean_raw = all_pmfs.mean(axis=0)
    ref_idx = int(np.argmin(mean_raw))
    ref_z = z_pos[ref_idx]

    anchored = all_pmfs - all_pmfs[:, ref_idx][:, None]   # anchor every replica at ref_z
    mean_pmf_ref = anchored.mean(axis=0)                  # pre-shift mean, referenced at ref_z
    se_pmf = anchored.std(axis=0, ddof=1) / np.sqrt(anchored.shape[0])  # se = 0 at ref_z

    # --- Single rigid shift so the MEAN curve is 0 in the bulk window ---
    # The shift is a constant, so it leaves the ref_z-anchored error bars (se_pmf)
    # unchanged; the same constant is applied to every batch curve so the
    # per-trajectory files stay on the same bulk = 0 scale as the summary.
    mask_win = (z_pos >= recal_lo) & (z_pos <= recal_hi)
    bulk_shift = mean_pmf_ref[mask_win].mean() if mask_win.any() else 0.0
    mean_pmf = mean_pmf_ref - bulk_shift
    print(f"{analog}: energy-minimum reference at z = {ref_z:.1f} A")

    # --- Write per-trajectory files: z | PMF_batch1 | PMF_batch2 | PMF_batch3 ---
    for traj in range(1, n_trajectories + 1):
        if traj not in traj_pmfs_map:
            continue
        shifted = [p - p[ref_idx] - bulk_shift for p in traj_pmfs_map[traj]]
        output_file = os.path.join(out_analog_dir, f"trajectory{traj}.dat")
        header = f"{'z':>12s} {'400-600ns':>12s} {'600-800ns':>12s} {'800-1000ns':>12s}"
        np.savetxt(
            output_file,
            np.column_stack([z_pos] + shifted),
            fmt="%12.6f",
            header=header,
            comments="",
        )
        print(f"Written {output_file}")

    # --- Summary: shifted mean and ref_z-anchored SE across all 9 batches ---
    summary_file = os.path.join(out_analog_dir, f"pmf_{analog}.dat")
    header = f"{'z':>12s} {'mean':>12s} {'se':>12s}"
    np.savetxt(
        summary_file,
        np.column_stack([z_pos, mean_pmf, se_pmf]),
        fmt="%12.6f",
        header=header,
        comments="",
    )
    print(f"Written {summary_file}")

    # --- Extra summary: mean + SE BEFORE the bulk shift (referenced at the
    #     energy-minimum depth ref_z, common to the 9 replicas) ---
    summary_ref_file = os.path.join(out_analog_dir, f"pmf_{analog}_emin.dat")
    np.savetxt(
        summary_ref_file,
        np.column_stack([z_pos, mean_pmf_ref, se_pmf]),
        fmt="%12.6f",
        header=header,
        comments="",
    )
    print(f"Written {summary_ref_file}")

print("Done.")


