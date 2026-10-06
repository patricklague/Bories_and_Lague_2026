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
- The two leaflets (`+z` and `−z`) are folded onto a single `|z|` grid.
- `C` is fixed by **recalibration**: `W` is shifted so that its mean over the
  bulk-water window `[recal_lo, recal_hi] = [40, 50] Å` is 0 (reference state =
  solute free in water).

The difficulty is that `p(z)` is estimated from a finite histogram, so bins that
were never visited have `p = 0` and `ln p = −∞`. The previous code avoided this
with `ln(p + ε)`, `ε = 1e-5`, which puts every empty bin at the **same** floor
`−k_B T ln ε`, independent of how much data was collected.

## 2. Bayesian (Dirichlet–multinomial) zero handling

Following Lovell, Chua & McGrath, *NAR Genomics and Bioinformatics* **2**(2),
lqaa040 (2020), the per-bin counts `n = (n_1, …, n_K)` are modelled as a
**multinomial** draw with a conjugate **Dirichlet(α)** prior. The posterior is
`Dirichlet(n_i + α)`, and its **posterior-mean proportion** is used in place of
the raw frequency:

```
p_i = (n_i + α) / (N + K·α)          N = Σ_i n_i ,  K = number of |z| bins
W_i = −k_B T · ln p_i  (then recalibrated)
```

Every `p_i > 0`, so the PMF is finite everywhere. An **empty** bin (`n_i = 0`)
no longer sits at a constant; it gets a **sample-size-aware floor**

```
p_floor = α / (N + K·α) ≈ α / N        W_floor ≈ −k_B T · ln(α / N)
```

i.e. "we saw 0 of N, so the probability is at most of order α/N." More sampling
(larger `N`) ⇒ deeper, more honest exclusion barriers, which is exactly the
count-scale dependence the paper argues must not be ignored.

Counts are recovered from the stored probability histogram using its own
quantisation: the smallest positive density equals one observation
(`quantum ≈ 1/N`), so `n_i = round(d_i / quantum)`. `N` can instead be given
explicitly through `N_samples`.

## 3. Parameters used and why

| Parameter | Value | Rationale |
|---|---|---|
| `k_B`, `T` | 0.008314 kJ/mol/K, 303.15 K | Simulation thermostat temperature; sets `k_B T ≈ 2.52 kJ/mol`. |
| `bayes_alpha` (α) | 0.5 | **Jeffreys prior** for the multinomial — the non-informative reference prior derived from Fisher information, invariant under reparametrisation. It adds only half a pseudo-count per bin, so it regularises empty bins while barely perturbing well-sampled ones. `α = 1` (Bayes–Laplace / uniform) is the main alternative; it biases slightly more toward the uniform composition. |
| `N_samples` | `None` (inferred) | With ~10⁵–10⁶ samples per batch the tails almost always contain single-count bins, so `1/min(positive density)` recovers `N` reliably. Set it explicitly if that assumption fails (see below). |
| `[recal_lo, recal_hi]` | `[40, 50] Å` | Bulk-water plateau far from the bilayer; defines `W = 0`. |
| `z_bin` | 1.0 Å | Grid resolution; also sets `K` and the counts-per-bin. |

With the current data (`N ≈ 5·10⁵`, `K ≈ 75`, `α = 0.5`): `K·α ≈ 37.5 ≪ N`, so
the prior is negligible for sampled bins (relative weight `K·α/N ≈ 7·10⁻⁵`) and
only acts where it is needed — the empty bins.

## 4. Criteria to check before adapting α (or N, bins)

The choice α = 0.5 is a principled default, **not** a tuned value. Verify the
following for your own system and change the parameters only if a check fails:

1. **Prior must stay negligible for sampled bins:** confirm `K·α ≪ N`
   (`N = Σ n_i` per batch). If `N` is small (short batches, few solutes, or a
   restrictive `monomer`/`multimer` filter), α begins to bias the well region —
   lower α or coarsen `z_bin` so bins hold more counts.
2. **α-sensitivity of the result:** recompute the PMF with α ∈ {0.5, 1, 1/K}.
   In the **sampled region** the curves should differ by far less than the
   block-to-block statistical error (`se` column). If they don't, the data are
   too sparse there and no prior choice can rescue it — collect more sampling.
3. **Physical plausibility of the floor:** `W_floor ≈ −k_B T ln(α/N)` should lie
   *above* the highest genuinely sampled barrier but not be absurdly large. If
   empty-bin walls dominate the figure, that region is simply unsampled; mark it
   rather than interpret it.
4. **Validity of the inferred N:** check that `1/min(positive density)` is close
   to an integer and matches the expected `frames_per_batch × n_solutes` (after
   any mode filter). If the smallest bin holds >1 count, inference
   underestimates `N`; pass the true value via `N_samples`.
5. **Recalibration window must be sampled:** the bins in `[recal_lo, recal_hi]`
   should contain real counts, otherwise the PMF zero is set by the floor rather
   than by the bulk-water reference. Widen the window or confirm bulk sampling.
6. **Bin-width convergence:** `z_bin` trades resolution against counts-per-bin.
   Check that the PMF in the region of interest is stable when `z_bin` is halved
   or doubled.
