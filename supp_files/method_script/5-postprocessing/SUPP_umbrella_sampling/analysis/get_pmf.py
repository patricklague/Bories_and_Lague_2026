#!/usr/bin/env python3
"""
get_pmf.py

Take the WHAM metadata produced by extract_metadata.py for three
independent 50 ns blocks of trajectory data (section101-section150,
section151-section200, section201-section250) and write, per block,
    <in-dir>/<block>/pmf_s1.dat        per-solute WHAM PMF (z, PMF, ...)
    <in-dir>/<block>/pmf_s2.dat
    <in-dir>/<block>/pmf-us-<analog>-<block>.dat
                                       block PMF (leaflets folded together)
and, combined over the three blocks,
    <out-dir>/pmf-us-<analog>.dat
        |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_blocks
    <out-dir>/pmf-us-<analog>-untranslated.dat
        same columns, block average before translation, SE = std/sqrt(n)

Protocol:
1) Within each block, solute 1 (+z leaflet, centers 0..+37 A) and solute 2
   (-z leaflet, centers -37..0 A) are run through WHAM separately; the
   solute-2 profile is mirrored to |z| and averaged with the solute-1 profile
   (ALC-2007 / MacCallum et al. two-leaflet treatment) without translation.
   This gives one PMF estimate per 50 ns block.
2) Each folded block is independently translated by subtracting its arithmetic
   mean over |z| >= --z-bulk (default 30 A), using bins sampled by all blocks.
3) The independently water-zeroed blocks are averaged.
4) Final SE comes from a bootstrap of the aligned blocks (exact enumeration
   of the n**n resamples when n**n <= --block-nboot), times sqrt(n/(n-1)).
   Referencing each block before resampling propagates its water-reference
   fluctuations and their covariance with the rest of the curve.
Each block has zero mean in the common water window, not necessarily zero
at every water bin. Arbitrary block constants cancel before averaging.

See get_pmf2.py for the alternative approach (single
50 ns block, error estimated from the two leaflets) used for shorter
(100 ns) trajectories.

By default the per-leaflet PMF is computed with binless WHAM (backend
"binless"): the window free energies are solved with the exact harmonic bias of
every sample (MBAR for a single unbiased state), and only the final unbiased
distribution is histogrammed. The output grid is 1 A (76 bins over -38..38 A,
edges on integer z), the 0.1 nm resolution at which MacCallum et al. (2008,
Biophys. J. 94:3393) report their PMFs (window spacing 0.1 nm, K = 3000
kJ mol^-1 nm^-2). Binned WHAM ("python", "grossfield") evaluates the bias at the
bin centre; with sigma = sqrt(kT/K) = 0.29 A, a 1 A bin then misplaces the bias
by several kT, and 0.25 A bins leave bin-scale noise. Override the output grid
with --num-bins or --bin-A.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# WHAM grid. Metadata contains window centers and force constants only; the
# PMF histogram grid is chosen here when WHAM is run.
# ---------------------------------------------------------------------------
DEFAULT_HIST_MIN = -38.0
DEFAULT_HIST_MAX =  38.0
DEFAULT_NUM_BINS = 76      # 1 A output bins (MacCallum et al. 0.1 nm resolution)
DEFAULT_TOL      = 1.0e-5
DEFAULT_TEMP     = 303.15
DEFAULT_NUM_PAD  = 0
DEFAULT_SEED     = 12345
DEFAULT_Z_BULK   = 30.0    # water reference window: |z| >= Z_BULK (to end of sampling)
DEFAULT_BLOCK_NBOOT = 10000
DEFAULT_BLOCK_SEED  = 2025
DEFAULT_MAX_ITER = 100000  # safety cap for the python WHAM iteration
WARM_START_BIN_A = 0.02    # fine binned WHAM used only to start the binless solver
BINLESS_TOL      = 1.0e-7  # binless fixed-point tolerance (reduced f); slow iteration needs it
KB_KCAL          = 0.0019872041   # Boltzmann constant in kcal/(mol K)


def _read_metadata(meta: Path) -> tuple[list[Path], np.ndarray, np.ndarray]:
    """Read Grossfield-style metadata: trajectory-file, window-center, K."""
    centers, kvals, files = [], [], []
    with open(meta) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            f, c, k = parts[0], float(parts[1]), float(parts[2])
            p = Path(f)
            if not p.is_absolute():
                p = (meta.parent / p).resolve()
            files.append(p)
            centers.append(c)
            kvals.append(k)
    return files, np.asarray(centers), np.asarray(kvals)


def wham_py(meta: Path, pmf: Path, log: Path,
            hist_min: float, hist_max: float, nbins: int,
            temp: float, tol: float, max_iter: int = DEFAULT_MAX_ITER) -> None:
    """Pure-numpy 1D WHAM. Assumes U_i(x) = (1/2) K_i (x - x0_i)^2
    (the colvars / NAMD convention; this is also the default in Grossfield
    WHAM >= 2.0.10).

    Writes a Grossfield-style PMF file:
        # coor  Free  +/-  Prob  +/-
        x_k     g_k   0    p_k   0
    and a short log file.
    """
    files, centers, kvals = _read_metadata(meta)
    nwin = len(files)
    kT = KB_KCAL * temp

    # The WHAM grid is independent of the umbrella window spacing. More bins
    # give a smoother/finer PMF output grid, but do not create new sampling.
    edges = np.linspace(hist_min, hist_max, nbins + 1)
    x = 0.5 * (edges[:-1] + edges[1:])                     # bin centres
    dx = (hist_max - hist_min) / nbins

    H = np.zeros((nwin, nbins), dtype=float)               # counts per window
    n = np.zeros(nwin, dtype=float)
    for i, f in enumerate(files):
        if not f.is_file() or f.stat().st_size == 0:
            continue
        data = np.loadtxt(f, comments=("#", "@"))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        vals = data[:, 1] if data.shape[1] >= 2 else data[:, 0]
        h, _ = np.histogram(vals, bins=edges)
        H[i] = h
        n[i] = h.sum()

    if nwin != 38:
        raise RuntimeError(f"expected 38 umbrella windows, found {nwin} in {meta}")
    if np.any(n == 0):
        empty = ", ".join(str(i) for i in np.where(n == 0)[0])
        raise RuntimeError(f"empty/out-of-range umbrella windows in {meta}: {empty}")
    if n.sum() == 0:
        raise RuntimeError(f"no samples landed in [{hist_min},{hist_max}]")

    normalized = H / n[:, None]
    adjacent_overlap = np.minimum(normalized[:-1], normalized[1:]).sum(axis=1)
    if np.any(adjacent_overlap == 0):
        pairs = ", ".join(
            f"{i}-{i + 1}" for i in np.where(adjacent_overlap == 0)[0]
        )
        raise RuntimeError(f"zero histogram overlap between windows {pairs} in {meta}")

    # bias matrix U_ik = (1/2) K_i (x_k - x0_i)^2 ;  W = exp(-U/kT)
    U = 0.5 * kvals[:, None] * (x[None, :] - centers[:, None]) ** 2
    W = np.exp(-U / kT)

    sum_H = H.sum(axis=0)                                  # total count per bin
    F = np.zeros(nwin)                                     # window free energies / kT
    rho = np.full(nbins, np.nan)

    occupied = sum_H > 0
    if not occupied.any():
        raise RuntimeError("all histogram bins are empty")

    with open(log, "w") as flog:
        flog.write(f"# python WHAM  T={temp}  bins={nbins}  "
                   f"range=[{hist_min},{hist_max}]  tol={tol}\n")
        flog.write(f"# {nwin} windows, total counts = {n.sum():.0f}\n")
        flog.write(f"# adjacent overlap: min={adjacent_overlap.min():.4f} "
               f"median={np.median(adjacent_overlap):.4f}\n")
        converged = False
        for it in range(1, max_iter + 1):
            # rho_k = sum_i H_ik  /  sum_i n_i exp(F_i) W_ik
            denom = (n * np.exp(F))[:, None] * W            # (nwin, nbins)
            denom_sum = denom.sum(axis=0)
            rho = np.where(occupied & (denom_sum > 0),
                           sum_H / np.where(denom_sum > 0, denom_sum, 1.0),
                           0.0)
            # F_i = -ln( sum_k W_ik rho_k dx )
            Z = (W * rho[None, :]).sum(axis=1) * dx
            with np.errstate(divide="ignore"):
                F_new = -np.log(np.where(Z > 0, Z, np.nan))
            F_new -= F_new[0]                              # gauge fix
            diff = np.nanmax(np.abs(F_new - F))
            F = F_new
            if it % 200 == 0 or diff < tol:
                flog.write(f"  iter {it:6d}  max|dF| = {diff:.3e}\n")
            if np.isfinite(diff) and diff < tol:
                flog.write(f"# converged in {it} iterations\n")
                converged = True
                break
        else:
            flog.write(f"# WARNING: not converged after {max_iter} iterations\n")

    if not converged:
        raise RuntimeError(f"WHAM did not converge for {meta} after {max_iter} iterations")

    g = np.where(rho > 0, -kT * np.log(np.where(rho > 0, rho, np.nan)), np.nan)
    if np.isfinite(g).any():
        g = g - np.nanmin(g)

    with open(pmf, "w") as fh:
        fh.write("#Coor\t\tFree\t+/-\t\tProb\t\t+/-\n")
        for xk, gk, pk in zip(x, g, rho):
            gs = f"{gk:.6f}" if np.isfinite(gk) else "nan"
            ps = f"{pk:.6e}" if np.isfinite(pk) and pk > 0 else "nan"
            fh.write(f"{xk:.6f}\t{gs}\t0\t{ps}\t0\n")


def wham_binless(meta: Path, pmf: Path, log: Path,
                 hist_min: float, hist_max: float, nbins: int,
                 temp: float, tol: float, max_iter: int = DEFAULT_MAX_ITER) -> None:
    """Binless 1D WHAM (MBAR for one unbiased state) with harmonic biases
    U_i(x) = (1/2) K_i (x - x0_i)^2.

    Window free energies f_i are solved self-consistently using the exact bias
    of every sample,  f_i = -ln sum_n exp(-u_i(x_n)) / sum_k N_k exp(f_k - u_k(x_n)),
    so no bias is evaluated at a bin centre. Each sample then carries the
    unbiased weight w_n = 1 / sum_k N_k exp(f_k - u_k(x_n)); the PMF is
    -kT ln of the weighted histogram on the output grid. Writes the same
    Grossfield-style file as wham_py.
    """
    files, centers, kvals = _read_metadata(meta)
    nwin = len(files)
    kT = KB_KCAL * temp
    if nwin != 38:
        raise RuntimeError(f"expected 38 umbrella windows, found {nwin} in {meta}")
    samples = []
    for f in files:
        if not f.is_file() or f.stat().st_size == 0:
            raise RuntimeError(f"missing/empty umbrella window file {f}")
        data = np.loadtxt(f, comments=("#", "@"))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        vals = data[:, 1] if data.shape[1] >= 2 else data[:, 0]
        if vals.size == 0 or not np.isfinite(vals).all():
            raise RuntimeError(f"empty or non-finite umbrella window {f}")
        samples.append(vals)
    n_k = np.array([len(s) for s in samples], dtype=float)
    x = np.concatenate(samples)

    edges = np.linspace(hist_min, hist_max, nbins + 1)
    xc = 0.5 * (edges[:-1] + edges[1:])

    # Overlap check and warm start from a fine binned WHAM (bias error << kT).
    fine_n = max(int(round((hist_max - hist_min) / WARM_START_BIN_A)), nbins)
    fine_edges = np.linspace(hist_min, hist_max, fine_n + 1)
    fine_x = 0.5 * (fine_edges[:-1] + fine_edges[1:])
    H = np.array([np.histogram(s, bins=fine_edges)[0] for s in samples], dtype=float)
    coarse = np.array([np.histogram(s, bins=edges)[0] for s in samples], dtype=float)
    if np.any(coarse.sum(axis=1) == 0):
        empty = ", ".join(str(i) for i in np.where(coarse.sum(axis=1) == 0)[0])
        raise RuntimeError(f"empty/out-of-range umbrella windows in {meta}: {empty}")
    normalized = H / H.sum(axis=1, keepdims=True)
    adjacent_overlap = np.minimum(normalized[:-1], normalized[1:]).sum(axis=1)
    if np.any(adjacent_overlap == 0):
        pairs = ", ".join(f"{i}-{i + 1}" for i in np.where(adjacent_overlap == 0)[0])
        raise RuntimeError(f"zero histogram overlap between windows {pairs} in {meta}")
    W_fine = np.exp(-0.5 * kvals[:, None] * (fine_x[None, :] - centers[:, None]) ** 2 / kT)
    n_fine = H.sum(axis=1)
    f = np.zeros(nwin)
    for _ in range(max_iter):
        rho = H.sum(axis=0) / np.maximum((n_fine * np.exp(f)) @ W_fine, 1e-300)
        f_new = -np.log(W_fine @ rho)
        f_new -= f_new[0]
        done = np.max(np.abs(f_new - f)) < tol
        f = f_new
        if done:
            break

    # Binless iterations; bias factors precomputed once (K x N).
    W = np.exp(-0.5 * kvals[:, None] * (x[None, :] - centers[:, None]) ** 2 / kT)
    converged = False
    with open(log, "w") as flog:
        flog.write(f"# binless WHAM (MBAR, single state)  T={temp}  output bins={nbins}  "
                   f"range=[{hist_min},{hist_max}]  tol={tol} (reduced units)\n")
        flog.write(f"# {nwin} windows, total samples = {x.size}\n")
        flog.write(f"# adjacent overlap ({WARM_START_BIN_A:g} A histogram): "
                   f"min={adjacent_overlap.min():.4f} median={np.median(adjacent_overlap):.4f}\n")
        flog.write(f"# warm start: binned WHAM with {fine_n} bins\n")
        for it in range(1, max_iter + 1):
            den = (n_k * np.exp(f)) @ W
            f_new = -np.log(W @ (1.0 / den))
            f_new -= f_new[0]
            diff = np.max(np.abs(f_new - f))
            f = f_new
            if it % 500 == 0 or diff < tol:
                flog.write(f"  iter {it:6d}  max|df| = {diff:.3e}\n")
            if np.isfinite(diff) and diff < tol:
                flog.write(f"# converged in {it} iterations\n")
                converged = True
                break
        else:
            flog.write(f"# WARNING: not converged after {max_iter} iterations\n")
    if not converged:
        raise RuntimeError(f"binless WHAM did not converge for {meta} after {max_iter} iterations")

    weights = 1.0 / ((n_k * np.exp(f)) @ W)
    p, _ = np.histogram(x, bins=edges, weights=weights)
    rho = p / np.diff(edges)
    g = np.where(rho > 0, -kT * np.log(np.where(rho > 0, rho, np.nan)), np.nan)
    if np.isfinite(g).any():
        g = g - np.nanmin(g)
    with open(pmf, "w") as fh:
        fh.write("#Coor\t\tFree\t+/-\t\tProb\t\t+/-\n")
        for xk, gk, pk in zip(xc, g, rho):
            gs = f"{gk:.6f}" if np.isfinite(gk) else "nan"
            ps = f"{pk:.6e}" if np.isfinite(pk) and pk > 0 else "nan"
            fh.write(f"{xk:.6f}\t{gs}\t0\t{ps}\t0\n")


def run_wham(wham_bin: str, hist_min: float, hist_max: float, nbins: int,
             tol: float, temp: float, num_pad: int, meta: Path,
             pmf: Path, log: Path, nboot: int, seed: int) -> None:
    """Run an external Grossfield WHAM binary with the same grid settings."""
    if shutil.which(wham_bin) is None:
        sys.exit(f"ERROR: '{wham_bin}' not on PATH. "
                 "Get it from http://membrane.urmc.rochester.edu/?page_id=126")
    cmd = [wham_bin, f"{hist_min}", f"{hist_max}", f"{nbins}", f"{tol}",
           f"{temp}", f"{num_pad}", str(meta), str(pmf)]
    if nboot > 0:
        cmd += [f"{nboot}", f"{seed}"]
    with open(log, "w") as fh:
        subprocess.run(cmd, cwd=meta.parent, stdout=fh, stderr=subprocess.STDOUT,
                       check=True)


def load_pmf(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load the coordinate and free-energy columns from a WHAM PMF file."""
    d = np.loadtxt(path, comments=("#", "@"))
    if d.ndim == 1:
        d = d.reshape(1, -1)
    return d[:, 0], d[:, 1]


def fold_side(zside: np.ndarray, gside: np.ndarray,
              target: np.ndarray) -> np.ndarray:
    """Interpolate one leaflet estimate onto the positive |z| output grid."""
    order = np.argsort(zside)
    zside, gside = zside[order], gside[order]
    finite = np.isfinite(gside)
    if finite.sum() < 2:
        return np.full_like(target, np.nan, dtype=float)
    return np.interp(target, zside[finite], gside[finite],
                     left=np.nan, right=np.nan)


def fold_leaflets(p1: Path, p2: Path
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fold the two offset solutes into one untranslated per-block PMF.

    In our setup solute 1 is restrained at 0..+37 A and solute 2 at -37..0 A.
    The two WHAM profiles are therefore independent estimates of the same
    |z|-PMF from opposite leaflets for one 50 ns block. Solute 2 is mirrored
    onto the positive grid and the two estimates are averaged WITHOUT any
    translation: each WHAM profile carries an arbitrary additive constant,
    and the average simply carries the mean of those constants, which is
    removed by the block's own water reference (see block_statistics).
    Bins not covered by both leaflets are set to NaN so that the average is
    never a mixture of curves with different constants.

    Returns |z| grid, block PMF, +z leaflet PMF, mirrored -z leaflet PMF.
    """
    z1, g1 = load_pmf(p1)
    z2, g2 = load_pmf(p2)
    if not np.allclose(z1, z2):
        sys.exit(f"{p1} and {p2} are on different z grids")
    z = z1

    pos = z >= 0
    neg = z < 0
    zp = z[pos]
    g_plus = g1[pos]
    g_minus = fold_side(-z[neg], g2[neg], zp)
    return zp, 0.5 * (g_plus + g_minus), g_plus, g_minus


def block_mean_se(pmfs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean and standard error over blocks (axis -2). pmfs: (..., n_blocks, n_z)."""
    n = pmfs.shape[-2]
    return pmfs.mean(axis=-2), pmfs.std(axis=-2, ddof=1) / np.sqrt(n)


def water_reference(mean: np.ndarray, se: np.ndarray,
                    mask_win: np.ndarray, se_floor: float = 1e-9) -> np.ndarray:
    """Inverse-variance weighted mean of the PMF over the water window.

    mean, se: (..., n_z). Returns a with shape (...). Falls back to an
    unweighted mean if any SE in the window is <= se_floor (degenerate
    bootstrap replicate in which every drawn block is the same block; such
    SEs are zero or floating-point noise).
    """
    m = mean[..., mask_win]
    s = se[..., mask_win]
    degenerate = np.any(s <= se_floor, axis=-1, keepdims=True)
    w = np.where(degenerate, 1.0, 1.0 / np.where(s > se_floor, s, 1.0) ** 2)
    return np.sum(w * m, axis=-1) / np.sum(w, axis=-1)


def block_statistics(z: np.ndarray, block_pmfs: np.ndarray, z_bulk: float,
                     n_boot: int, seed: int, small_sample_correction: bool = True
                     ) -> dict:
    """Independently water-zero blocks, then average and bootstrap them.

    Every block uses the same water bins, sampled by all blocks. The raw
    mean/SE are retained solely for the untranslated diagnostic output.
    When n_blocks**n_blocks <= n_boot, all equally likely resamples are
    enumerated (exact bootstrap) instead of drawing them at random.
    """
    n_blocks = block_pmfs.shape[0]
    if n_blocks < 2:
        sys.exit("ERROR: at least 2 blocks are needed to estimate the SE")

    valid = np.all(np.isfinite(block_pmfs), axis=0)
    mask_win = valid & (z >= z_bulk)
    if not mask_win.any():
        sys.exit(f"ERROR: no bin sampled by every block in |z| >= {z_bulk} A")
    P = block_pmfs[:, valid]
    win = mask_win[valid]

    mean_raw, se_raw = block_mean_se(P)
    block_references = P[:, win].mean(axis=1)
    aligned = P - block_references[:, None]
    mean = aligned.mean(axis=0)
    a = float(block_references.mean())

    if n_blocks ** n_blocks <= n_boot:
        idx = np.array(np.meshgrid(*[np.arange(n_blocks)] * n_blocks,
                                   indexing="ij")).reshape(n_blocks, -1).T
        method = f"exact (all {len(idx)} resamples)"
    else:
        idx = np.random.default_rng(seed).integers(0, n_blocks, size=(n_boot, n_blocks))
        method = f"{n_boot} random resamples, seed {seed}"
    boot_shifted = aligned[idx].mean(axis=1)
    boot_a = block_references[idx].mean(axis=1)

    corr = np.sqrt(n_blocks / (n_blocks - 1)) if small_sample_correction else 1.0
    # Population std over the resamples: with exact enumeration this is the
    # exact bootstrap SE; with random resamples the difference is negligible.
    se_final = corr * boot_shifted.std(axis=0)
    sigma_a = float(corr * boot_a.std())

    def on_grid(v: np.ndarray) -> np.ndarray:
        out = np.full(z.shape, np.nan)
        out[valid] = v
        return out

    return {
        "a": a, "sigma_a": sigma_a, "block_references": block_references,
        "mean_raw": on_grid(mean_raw), "se_raw": on_grid(se_raw),
        "mean": on_grid(mean), "se": on_grid(se_final),
        "n_blocks": n_blocks, "n_win": int(mask_win.sum()), "method": method,
        "corr": corr,
    }


def write_block(out: Path, z: np.ndarray, g: np.ndarray, g_plus: np.ndarray,
                g_minus: np.ndarray, a: float, analog: str, block: str,
                z_bulk: float, bin_width: float) -> None:
    """Write one block PMF shifted by its own water reference a.

    The leaflet SE column is a diagnostic only (leaflets each zeroed on the
    water window before comparison, since their WHAM constants are arbitrary);
    it is not propagated to the final uncertainty.
    """
    win = (z >= z_bulk) & np.isfinite(g_plus) & np.isfinite(g_minus)
    with np.errstate(invalid="ignore"):
        dp = g_plus - (np.mean(g_plus[win]) if win.any() else np.nan)
        dm = g_minus - (np.mean(g_minus[win]) if win.any() else np.nan)
        se = np.abs(dp - dm) / 2.0          # = std(ddof=1)/sqrt(2) for 2 values
    n = np.isfinite(g_plus).astype(int) + np.isfinite(g_minus).astype(int)
    with open(out, "w") as fh:
        fh.write(f"# PMF for analog '{analog}', block '{block}' "
                 f"(ALC-2007 two-leaflet averaging, {bin_width:g} A bins)\n")
        fh.write("# columns: |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_leaflets\n")
        fh.write(f"# folded block shifted by its own arithmetic-mean water "
                 f"reference a = {a:.5f} kcal/mol (|z| >= {z_bulk:g} A, "
                 "bins sampled by all blocks)\n")
        fh.write("# SE = leaflet-to-leaflet spread, diagnostic only -- "
                 "final uncertainty comes from the block bootstrap\n")
        for zi, gi, si, ni in zip(z, g - a, se, n):
            fh.write(f"{zi:8.3f} {gi:12.5f} {si:12.5f} {int(ni):4d}\n")


def write_combined(out: Path, z: np.ndarray, mean: np.ndarray, se: np.ndarray,
                   n_blocks: int, block_names: list[str], analog: str,
                   header_lines: list[str]) -> int:
    with open(out, "w") as fh:
        fh.write(f"# PMF for analog '{analog}'  "
                 f"(averaged over {len(block_names)} 50 ns blocks: "
                 f"{', '.join(block_names)})\n")
        fh.write("# columns: |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_blocks\n")
        for line in header_lines:
            fh.write(f"# {line}\n")
        for zi, gi, si in zip(z, mean, se):
            fh.write(f"{zi:8.3f} {gi:12.5f} {si:12.5f} {n_blocks:4d}\n")
    return int(np.isfinite(mean).sum())

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analog",   required=True)
    ap.add_argument("--in-dir",   required=True, type=Path,
                    help="parent directory containing one subdirectory per "
                         "50 ns block (see --blocks), each with "
                         "metadata_s1.dat/metadata_s2.dat from extract_metadata.py")
    ap.add_argument("--out-dir",  type=Path, default=None,
                    help="where to write the combined pmf-us-<analog>.dat "
                         "(default: same as --in-dir)")
    ap.add_argument("--blocks", nargs="+",
                    default=["block1", "block2", "block3"],
                    help="block subdirectory names under --in-dir, in order "
                         "(default: block1 block2 block3, matching "
                         "section101-150 / section151-200 / section201-250)")
    ap.add_argument("--backend", choices=["binless", "python", "grossfield"],
                    default="binless",
                    help="WHAM implementation (default: binless, exact per-sample "
                         "bias; python/grossfield are binned WHAM)")
    ap.add_argument("--wham-bin", default="wham",
                    help="path/name of Grossfield wham binary (--backend grossfield)")
    ap.add_argument("--hist-min", type=float, default=DEFAULT_HIST_MIN)
    ap.add_argument("--hist-max", type=float, default=DEFAULT_HIST_MAX)
    ap.add_argument("--num-bins", type=int, default=DEFAULT_NUM_BINS,
                    help=f"number of WHAM histogram bins (default {DEFAULT_NUM_BINS})")
    ap.add_argument("--bin-A",    type=float, default=None,
                    help="bin width in A; overrides --num-bins when set")
    ap.add_argument("--tol",      type=float, default=DEFAULT_TOL)
    ap.add_argument("--temp",     type=float, default=DEFAULT_TEMP)
    ap.add_argument("--num-pad",  type=int,   default=DEFAULT_NUM_PAD)
    ap.add_argument("--z-bulk",   type=float, default=DEFAULT_Z_BULK,
                    help="water reference window is |z| >= z-bulk up to the "
                         "end of sampling (default 30 A)")
    ap.add_argument("--nboot",    type=int,   default=0,
                    help="Grossfield WHAM internal bootstrap rounds (not used "
                         "for the final SE)")
    ap.add_argument("--seed",     type=int,   default=DEFAULT_SEED)
    ap.add_argument("--block-nboot", type=int, default=DEFAULT_BLOCK_NBOOT,
                    help="block-bootstrap replicates for the final SE; all "
                         "n_blocks**n_blocks resamples are enumerated when "
                         "that is smaller (default 10000)")
    ap.add_argument("--block-seed",  type=int, default=DEFAULT_BLOCK_SEED)
    args = ap.parse_args()

    out_dir = args.out_dir or args.in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    hist_range = args.hist_max - args.hist_min
    if hist_range <= 0:
        sys.exit("ERROR: --hist-max must be greater than --hist-min")
    if args.bin_A is not None:
        if args.bin_A <= 0:
            sys.exit("ERROR: --bin-A must be positive")
        nbins = int(round(hist_range / args.bin_A))
        bin_width = hist_range / nbins
    else:
        if args.num_bins <= 0:
            sys.exit("ERROR: --num-bins must be positive")
        nbins = args.num_bins
        bin_width = hist_range / nbins

    block_z: np.ndarray | None = None
    block_means: list[np.ndarray] = []
    block_leaflets: list[tuple[np.ndarray, np.ndarray]] = []
    used_blocks: list[str] = []
    for block in args.blocks:
        block_dir = args.in_dir / block
        meta1 = block_dir / "metadata_s1.dat"
        meta2 = block_dir / "metadata_s2.dat"
        missing = [m for m in (meta1, meta2) if not m.is_file() or m.stat().st_size == 0]
        if missing:
            print(f"  ! block '{block}': missing/empty metadata "
                  f"({', '.join(str(m) for m in missing)}), skipped")
            continue

        pmf1 = block_dir / "pmf_s1.dat"
        pmf2 = block_dir / "pmf_s2.dat"
        log1 = block_dir / "wham_s1.log"
        log2 = block_dir / "wham_s2.log"

        print(f"  [{block}] WHAM s1: bins={nbins} ({bin_width:g} A) "
              f"range=[{args.hist_min},{args.hist_max}]  T={args.temp}  "
              f"backend={args.backend}")
        if args.backend == "binless":
            tol = min(args.tol, BINLESS_TOL)
            wham_binless(meta1, pmf1, log1, args.hist_min, args.hist_max, nbins,
                         args.temp, tol)
            print(f"  [{block}] WHAM s2: bins={nbins}")
            wham_binless(meta2, pmf2, log2, args.hist_min, args.hist_max, nbins,
                         args.temp, tol)
        elif args.backend == "python":
            wham_py(meta1, pmf1, log1, args.hist_min, args.hist_max, nbins,
                    args.temp, args.tol)
            print(f"  [{block}] WHAM s2: bins={nbins}")
            wham_py(meta2, pmf2, log2, args.hist_min, args.hist_max, nbins,
                    args.temp, args.tol)
        else:
            run_wham(args.wham_bin, args.hist_min, args.hist_max, nbins, args.tol,
                     args.temp, args.num_pad, meta1, pmf1, log1, args.nboot, args.seed)
            print(f"  [{block}] WHAM s2: bins={nbins}")
            run_wham(args.wham_bin, args.hist_min, args.hist_max, nbins, args.tol,
                     args.temp, args.num_pad, meta2, pmf2, log2, args.nboot, args.seed)

        # Step 1: untranslated per-block PMF (leaflets folded together)
        zp, mean, g_plus, g_minus = fold_leaflets(pmf1, pmf2)
        if block_z is None:
            block_z = zp
        elif not np.allclose(block_z, zp):
            sys.exit(f"ERROR: block '{block}' PMF grid differs from earlier blocks")
        block_means.append(mean)
        block_leaflets.append((g_plus, g_minus))
        used_blocks.append(block)

    if not block_means:
        sys.exit("ERROR: no usable blocks found under " + str(args.in_dir))

    # Steps 2-4: independently water-zero blocks, average, and bootstrap.
    st = block_statistics(block_z, np.vstack(block_means), args.z_bulk,
                          args.block_nboot, args.block_seed)
    print(f"  mean of block water references = {st['a']:.5f} +/- "
          f"{st['sigma_a']:.5f} kcal/mol (bootstrap), "
          f"{st['n_win']} common bins in |z| >= {args.z_bulk:g} A")

    # Step 5: outputs
    for block, mean, (g_plus, g_minus), a in zip(
            used_blocks, block_means, block_leaflets, st["block_references"]):
        block_pmf = args.in_dir / block / f"pmf-us-{args.analog}-{block}.dat"
        write_block(block_pmf, block_z, mean, g_plus, g_minus, float(a),
                    args.analog, block, args.z_bulk, bin_width)
        print(f"  [{block}] -> {block_pmf}")

    ref_lines = [
        f"each folded block independently shifted by its arithmetic-mean PMF "
        f"over |z| >= {args.z_bulk:g} A ({st['n_win']} bins sampled by all blocks)",
        f"mean of block references = {st['a']:.5f} +/- {st['sigma_a']:.5f} "
        "kcal/mol (bootstrap; references depend on arbitrary WHAM constants)",
    ]
    combined = out_dir / f"pmf-us-{args.analog}.dat"
    nbin_ok = write_combined(
        combined, block_z, st["mean"], st["se"], st["n_blocks"], used_blocks,
        args.analog,
        ref_lines + [
            f"per-leaflet WHAM backend = {args.backend}, output bin width = {bin_width:g} A",
            "PMF = average of independently water-zeroed blocks; each block folds solute 1 (+z) and "
            "mirrored solute 2 (-z) leaflets",
            f"SE = block bootstrap of independently water-zeroed blocks ({st['method']}), "
            f"x sqrt(n/(n-1)) = {st['corr']:.4f} small-sample correction",
        ])
    print(f"==> {combined}   ({nbin_ok} usable bins, {len(used_blocks)} blocks)")

    untranslated = out_dir / f"pmf-us-{args.analog}-untranslated.dat"
    write_combined(
        untranslated, block_z, st["mean_raw"], st["se_raw"], st["n_blocks"],
        used_blocks, args.analog,
        ["untranslated block average BEFORE independent block water referencing "
         "(arbitrary WHAM constants); SE = std/sqrt(n) across raw blocks"])
    print(f"    {untranslated}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
