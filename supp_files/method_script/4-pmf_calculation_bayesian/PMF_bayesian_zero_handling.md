# PMF from density histograms with Bayesian zero handling

This note documents the equations and parameters used in
[`pmf-from-distribution.py`](./pmf-from-distribution.py).

## 1. PMF by Boltzmann inversion

For an equilibrium ensemble the probability of finding the solute COM at depth
`z` is `p(z) ∝ exp(−W(z)/k_B T)`, so the potential of mean force (PMF) is

```
W(z) = −k_B T · ln p(z) + C
```

- `k_B = 0.008314 kJ mol⁻¹ K⁻¹`, `T = 303.15 K`  →  `k_B T ≈ 2.52 kJ mol⁻¹`.
- The two leaflets (`+z` and `−z`) are folded onto a single `|z|` grid by summing
  their counts.
- `C` is set by the **water reference** `a` (section 4): `W` is shifted so that
  the bulk-water window `[recal_lo, recal_hi] = [40, 50] Å` is at 0 (reference
  state = solute free in water).

The difficulty is that `p(z)` is estimated from a finite histogram, so bins that
were never visited have `p = 0` and `ln p = −∞`. The previous code avoided this
with `ln(p + ε)`, `ε = 1e-5`, which puts every empty bin at the **same** floor
`−k_B T ln ε`, independent of how much data was collected. With `N = 520 000`,
`ε` corresponds to ~5 counts, i.e. an empty bin was treated as *more* probable
than a bin observed once.

## 2. Bayesian (Dirichlet–multinomial) zero handling

Following Lovell, Chua & McGrath, *NAR Genomics and Bioinformatics* **2**(2),
lqaa040 (2020), the per-bin counts `n = (n_1, …, n_K)` are modelled as a
**multinomial** draw with a conjugate **Dirichlet(α)** prior. The posterior is
`Dirichlet(n_i + α)`, and its **posterior-mean proportion** is used in place of
the raw frequency:

```
p_i = (n_i + α) / (N + K·α)          N = Σ_i n_i ,  K = number of |z| bins
W_i = −k_B T · ln p_i  (absolute; no translation at this stage)
```

Every `p_i > 0`, so the PMF is finite everywhere. An **empty** bin (`n_i = 0`)
no longer sits at a constant; it gets a **sample-size-aware floor**

```
p_floor = α / (N + K·α) ≈ α / N        W_floor ≈ −k_B T · ln(α / N)
```

i.e. "we saw 0 of N, so the probability is at most of order α/N." More sampling
(larger `N`) ⇒ deeper, more honest exclusion barriers, which is exactly the
count-scale dependence the paper argues must not be ignored. The transform is
applied to **each 200 ns block separately**, using only that block's counts, so
the nine block PMFs remain statistically independent.

**Counts.** The stored histograms are normalised probabilities printed with 6
decimals. They are converted back to exact integer counts, `n_i = round(d_i·N)`,
with `N` = number of molecule-frames in the block, taken from the raw data
(`raw_data_summary.dat`, otherwise the line count of the raw file; 520 000 per
block, 260 000 for the `-1` analogs). The rounded counts then sum exactly to `N`
(residual ≤ 0.24 count). The earlier rule `N = 1/min(positive density)` is
**not** used: when no bin holds a single count it underestimated `N` by up to a
factor 6. A warning is printed if the residual exceeds 0.45 count or the total
differs from `N` by more than 1 % (currently only `scrn-1`, trajectory 2,
block 3, whose stored distribution differs slightly from the raw data).

## 3. Block averaging, water reference and uncertainties

Per analog, with 3 trajectories × 3 blocks of 200 ns (`n = 9` blocks):

1. **Block PMFs** `W_b(z)` from section 2, without any translation.
2. **Block average** `W̄(z) = mean_b W_b(z)` and `SE(z) = std_b(W_b, ddof=1)/√n`,
   still without translation.
3. **Water reference** `a` = inverse-variance weighted mean of `W̄` over
   `[40, 50] Å`, with weights `1/SE²`. If an SE in the window is ≤ 1e-9 kJ/mol
   (identical values in every block, e.g. a window empty everywhere), an
   unweighted mean is used.
4. **Translation** `W'(z) = W̄(z) − a`. Its uncertainty comes from a **block
   bootstrap**: the nine block PMFs are resampled with replacement (10 000
   replicates), and steps 2–4 are repeated on each replicate. The SE of `W'` is
   the standard deviation of the shifted replicates × `√(n/(n−1))`. This
   propagates the uncertainty of `a` with its full covariance with `W̄`. The
   error is therefore smaller inside the reference window, where `W̄` and `a`
   are correlated, and includes `σ_a` elsewhere. A simple quadrature sum
   `√(SE² + σ_a²)` would be wrong, because the probability normalisation makes
   the whole water level of a block move together (the bootstrap `σ_a` is
   ~3× the independent-bin value `1/√Σw`).
5. **Outputs** in `total/<analog>/`:
   - `pmf_<analog>.dat`: `z`, `mean` (= `W'`), `se`, `n_empty` (number of blocks
     with zero counts in the bin; `n_empty = 9` means the value is the prior
     floor only);
   - `pmf_<analog>_untranslated.dat`: `z`, `W̄`, `SE` (step 2);
   - `trajectory<n>.dat`: block PMFs shifted by the common `a`;
   - `total/water_reference.dat`: `a`, bootstrap `σ_a` and the
     independent-bin estimate for every analog.

   The former `pmf_<analog>_emin.dat` (reference at the energy minimum) is no
   longer produced.

## 4. Parameters used and why

| Parameter | Value | Rationale |
|---|---|---|
| `k_B`, `T` | 0.008314 kJ/mol/K, 303.15 K | Simulation thermostat temperature; sets `k_B T ≈ 2.52 kJ/mol`. |
| `bayes_alpha` (α) | 0.5 | **Jeffreys prior** for the multinomial: the non-informative reference prior derived from Fisher information, invariant under reparametrisation. It adds only half a pseudo-count per bin, so it regularises empty bins while barely perturbing well-sampled ones. `α = 1` (Bayes–Laplace / uniform) is the main alternative; it biases slightly more toward the uniform composition. |
| `N_samples` | `None` (from raw data) | Per-block sample size. `None` reads it from the raw data (section 2); set a value only to override it. |
| `n_trajectories`, `n_batches` | 3, 3 | Three independent trajectories, each split into three non-overlapping 200 ns blocks (400–1000 ns), treated as independent samples (assumed; the correlation time of the solute depth was not measured). |
| `[recal_lo, recal_hi]` | `[40, 50] Å` | Bulk-water plateau far from the bilayer; defines `W' = 0`. |
| `z_bin` | 1.0 Å | Grid resolution; also sets `K` and the counts per bin. |
| `n_bootstrap` | 10 000 | Monte Carlo noise of the bootstrap SE < 1 % (0.7 % spread across seeds; 2.4 % at 1 000 and 4.7 % at 200 replicates). |
| `bootstrap_seed` | 2025 | Arbitrary; only makes the output reproducible (another seed changes the SE by ~0.7 %). |
| `small_sample_correction` | `True` | Multiplies the bootstrap SE by `√(n/(n−1)) = 1.061` (n = 9). This removes the known downward bias of the bootstrap variance of a mean, so the SE matches the conventional `s/√n`; the bootstrap/analytical SE ratio is 0.99–1.03. |

With the current data (`N = 5.2·10⁵`, `K ≈ 75`, `α = 0.5`): `K·α ≈ 37.5 ≪ N`, so
the prior is negligible for sampled bins (relative weight `K·α/N ≈ 7·10⁻⁵`) and
only acts where it is needed: the empty bins.

## 5. Criteria to check before adapting α (or N, bins)

The choice α = 0.5 is a principled default, **not** a tuned value. Verify the
following for your own system and change the parameters only if a check fails:

1. **Prior must stay negligible for sampled bins:** confirm `K·α ≪ N`
   (`N = Σ n_i` per batch). If `N` is small (short batches, few solutes, or a
   restrictive `monomer`/`multimer` filter), α begins to bias the well region:
   lower α or coarsen `z_bin` so bins hold more counts.
2. **α-sensitivity of the result:** recompute the PMF with α ∈ {0.5, 1, 1/K}.
   In the **sampled region** the curves should differ by far less than the
   block-to-block statistical error (`se` column). If they don't, the data are
   too sparse there and no prior choice can rescue it: collect more sampling.
   *Current data:* the final (shifted) PMF changes by ≤ 0.3 kJ/mol for polar
   and charged analogs (sca, scs, scd), mostly in sparse tail bins. For
   hydrophobic analogs whose water window is partly empty, it changes by
   1.0–1.6 kJ/mol (α = 1) and 5–8 kJ/mol (α = 1/K) for scl and scm-1. Their
   water reference therefore depends on the prior, and this check fails
   for them.
3. **Physical plausibility of the floor:** `W_floor ≈ −k_B T ln(α/N)` should lie
   *above* the highest genuinely sampled barrier but not be absurdly large. If
   empty-bin walls dominate the figure, that region is simply unsampled; mark it
   (`n_empty` column) rather than interpret it.
4. **Validity of N:** the script recovers exact counts from the raw-data sample
   size and warns if they are not integers or do not sum to `N`. For filtered
   modes (`monomer`/`multimer`) the per-block totals differ from the raw line
   count, and this reconstruction does not apply.
5. **Recalibration window must be sampled:** the bins in `[recal_lo, recal_hi]`
   should contain real counts, otherwise the PMF zero is set by the floor rather
   than by the bulk-water reference. The script warns when cells of the window
   are empty. Currently this concerns 12 hydrophobic analogs. `sci-1` has no
   count at all in the window, so its reference is pure prior and its
   `σ_a = 0` is not a real uncertainty.
6. **Bin-width convergence:** `z_bin` trades resolution against counts-per-bin.
   Check that the PMF in the region of interest is stable when `z_bin` is halved
   or doubled.