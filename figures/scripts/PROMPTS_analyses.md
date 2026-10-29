# Analysis prompts for Claude — outstanding items from TODO2.md

**How to use this file.** Each numbered section below is a self-contained prompt you can paste into a fresh Claude conversation. Fill in every `[[…]]` placeholder with the actual paths, module names, or numeric values from your analysis repo before sending. Read the shared block *Conventions and inputs* first — every analysis prompt assumes it.

---

## Conventions and inputs (shared, must precede every prompt)

Copy this block once at the top of the conversation, or include it verbatim inside every individual prompt you send.

> **Study.** MD survey of amino acid side-chain analog insertion into a POPC bilayer. Submission target: BBA Biomembranes.
>
> **Coordinate convention.**
> - `z` is the bilayer normal, with `z = 0` at the instantaneous bilayer centre of mass. Trajectories were recentred prior to analysis.
> - Analog position is the analog's centre of mass.
> - Bilayer regions (operational, POPC-specific): Region I `|z| < 9.5 Å`; Region II `9.5 ≤ |z| < 20 Å`; Region III `20 ≤ |z| ≤ 30 Å`; Region IV `|z| > 30 Å`.
> - Where properties are averaged over both leaflets, `z` is folded to `|z|`. Where leaflets are analysed separately (e.g. leaflet-symmetry PMF), the fold is **disabled** and both signs of `z` are retained.
>
> **PMF convention.**
> - Histogram bin width `Δz = 1.0 Å`. Density `ρ(z)` obtained by converting the histogram counts into a normalised probability density.
> - `T = 303.15 K`. `k_B T ≈ 2.522 kJ/mol`.
> - `PMF(z) = -k_B T · ln ρ(z) + C`, with `C` chosen so `⟨PMF⟩ = 0` over the bulk-water region `40 ≤ |z| ≤ 50 Å`.
> - Bins with 0 samples are flagged "sampling-limited" and rendered as dashed segments in plots and excluded from downstream summary statistics.
>
> **Mean and standard-error convention (unbiased-MD analyses).** Consistent with the Methods section of the current manuscript:
> - Production window: **final 600 ns** of each of **3 independent trajectories** per system.
> - Block size: **200 ns**, non-overlapping. This yields **9 blocks total per system** (3 blocks × 3 replicas).
> - Block value `x_i` is the analysis output computed on block `i` in isolation (not a running mean across blocks).
> - Reported **mean** `x̄ = (1/N) Σ x_i` with `N = 9`.
> - Reported **standard error** `SE = s / √N` with `s = np.std(x, ddof=1)` (sample standard deviation with N−1 in the denominator).
> - For **paired** analyses (POPC-only vs analog-containing systems), compute nine paired **block differences** `Δ_i = x_i^{analog} − x_i^{POPC-only}` using matched block indices across replicas, then take mean and SE of `{Δ_i}` with the same formula.
>
> **Mean and standard-error convention (US analyses).**
> - Production window: **final 150 ns** of each 250-ns US window (first 100 ns discarded as additional equilibration).
> - Block size: **50 ns** → **3 blocks per window**.
> - PMF SE from WHAM: run WHAM on each block independently, then compute `SE(z)` across the 3 blocks with the same `s / √N` formula (`N = 3`).

---

## Analysis A1 — Leaflet-symmetry PMF as a convergence probe

**TODO2 link:** §1.1 (per-analog leaflet-symmetry metric).

> **Objective.** For every analog, compute the PMF separately from the upper and lower leaflets of the POPC bilayer and quantify the leaflet asymmetry as an independent convergence probe. A converged PMF must be symmetric about `z = 0` within its statistical uncertainty.
>
> **Inputs.**
> - Density data files : `figures/data/distribution_data/total/{resid}/` - 1 file per trajectory (separate in 200ns)
> - Analog list: `[[ANALOG LIST: ALA, VAL, LEU, ILE, MET, CYS0, SER, THR, ASN, GLN, PHE, TYR0, TRP, PRO, GLYD, ASP0, GLU0, HSE0, HSD0, LYS0, ARG0, ASP-, GLU-, CYS-, TYR-, HSP+, LYS+, ARG+]]`
>
> **Method.**
> 1. Boltzmann-invert each histogram to `PMF_up(z)` and `PMF_down(z)` using the shared PMF convention. Reference each independently to zero over its own bulk-water region (`z ∈ [40, 50] Å` for `PMF_up`, `z ∈ [-50, -40] Å` for `PMF_down`).
> 2. Define the leaflet-asymmetry residual `ΔG_sym(z) = PMF_up(z) − PMF_down(−z)` on the common depth grid `z ∈ [0, 50] Å`, evaluated only in bins where both branches pass the `BIN_COUNT_MIN` threshold.
> 3. Compute per-analog scalar summaries:
>    - `max|ΔG_sym|` — the largest absolute asymmetry within the well-sampled depth range,
>    - `RMS(ΔG_sym)` — the root-mean-square asymmetry across well-sampled bins,
>    - the depth `z*` at which `max|ΔG_sym|` occurs.
> 4. Estimate the mean and SE of these three scalars using the standard **9-block convention**: compute `ΔG_sym(z)`, `max|ΔG_sym|`, `RMS(ΔG_sym)`, and `z*` independently on each of the 9 200-ns blocks; then take mean and SE across the 9 blocks with `s / √N`, `N = 9`.
> 5. Classify each analog: **pass** if `max|ΔG_sym|` mean + 2·SE is below the tolerance `2 k_B T ≈ 5 kJ/mol`, **flag** otherwise.
>
> **Outputs.**
> - `rev2/leaflet_symmetry/leaflet_symmetry_summary.csv` with one row per analog: `analog, max_abs_dGsym_mean, max_abs_dGsym_se, rms_dGsym_mean, rms_dGsym_se, zstar_mean, zstar_se, pass/flag`.
> - `rev2/leaflet_symmetry/per_analog_dGsym_curves.csv` with one row per (analog, z) grid point and columns `dGsym_mean, dGsym_se`.
> - `rev2/leaflet_symmetry/figures/dGsym_vs_z_<ANALOG>.pdf` — one figure per analog showing `PMF_up(z)`, `PMF_down(−z)`, and `ΔG_sym(z)` with ±SE bands, plus dashed segments for sampling-limited bins.
> - A short markdown report `rev2/leaflet_symmetry/REPORT.md` summarising the pass/flag classification and naming the analogs that fail the tolerance.

---

## Analysis A2 — Cumulative running-mean depth-distribution convergence

**TODO2 link:** §1.1 (redraw supplemental convergence panels; y-axis = cumulative running mean, not instantaneous population).

> **Objective.** For each analog, redo the current "populations vs time" plots (Fig. S7 successor) with a **cumulative running-mean** representation, so a converged trajectory settles inside its ±SE band and Reviewer 2's "still drifting" complaint can be visually rebutted or accepted.
>
> **Inputs.**.
> - COM data files : `figures/data/distribution_data/raw_data/{resid}/*` - 1 file per trajectory
> - Analog list: `[[ANALOG LIST: LEU, CYS, SER, TYR, ASP-,ARG+]]`

> **Method.**
> 1. For each of the 3 replicas and each analog, compute the per-frame Region-I…IV occupancy indicator `r_k(t)`, `k ∈ {I, II, III, IV}`, defined as the number of analog molecules whose folded `|z(t)|` falls in region `k` at frame `t`. Save at the trajectory's frame stride (10 ps).
> 2. For each replica and region, form the **cumulative running mean** `R_k(t) = (1/t) ∫_0^t r_k(t') dt'`. Return this on a coarser output grid at 1 ns resolution.
> 3. Estimate the running SE using a **block-based estimator**: at each output time `t ≥ 200 ns`, use non-overlapping 200-ns blocks of `r_k` from `[0, t]`, take the mean of each block, and compute `SE(t) = s(blocks) / √N_blocks(t)` where `N_blocks(t) = floor(t / 200 ns)`. Values with `N_blocks(t) < 2` are undefined.
> 4. For the final-window value that appears in the manuscript, use the **standard 9-block convention** on the final 600 ns only, unchanged from the manuscript methodology.
>
> **Outputs.**
> - `rev2/running_mean/figures/running_mean_<ANALOG>.png` — one figure per analog with 4 panels (one per region), each showing 3 replica curves + ±SE band from step 3, a horizontal line at the final-600-ns mean, and a vertical dashed line at 400 ns marking the equilibration/production boundary.
> - `rev2/running_mean/final_600ns_values.csv` — one row per analog with `region_I_mean, region_I_se, …, region_IV_se` on the 9-block convention.
> - `rev2/running_mean/REPORT.md` — flag any (analog, region) pair whose final-600 ns mean lies outside the running ±SE band of the last 100 ns of any replica; those are the ones Reviewer 2 would call "still drifting".

---

## Analysis A3 — Region-population running-average vs block index (companion to block-size dependence)

**TODO2 link:** §1.1 last bullet.

> **Objective.** Answer the "oscillation around the final value" request literally. For a representative set of analogs, plot the block-average Region-I…IV populations against block index, so oscillation around a plateau is directly visible.
>
> **Inputs.** Same as A2. Analog subset: `[[ANALOG SUBSET FOR THIS PLOT: e.g. LEU, CYS0, SER, TYR0, ARG+, ASP-]]` (matches the current Fig. S7 subset).
>
> **Method.**
> 1. Use the same 9-block partitioning as the manuscript.
> 2. Per analog and per region, compute the block-mean population `x_i^{(analog, region)}` for `i ∈ 1…9`.
> 3. Plot `x_i` versus block index `i`, with different markers or colours for the three replicas (blocks 1–3 from replica 1, 4–6 from replica 2, 7–9 from replica 3).
> 4. Overlay a horizontal band `x̄ ± SE` computed on the same 9-block convention.
>
> **Outputs.**
> - `rev2/block_index/figures/block_index_<ANALOG>.png` — 4 sub-panels per analog (regions I–IV).
> - `rev2/block_index/oscillation_amplitude.csv` — one row per (analog, region): `mean, se, block_range = max_i(x_i) − min_i(x_i)`. Analogs whose `block_range > 2·SE` in any region are flagged.

---

## Analysis A4 — US window-overlap and WHAM decorrelation diagnostics

**TODO2 link:** §1.4 (respond to Reviewer 2's question on the US force constant).

> **Objective.** Report the observed window-to-window histogram overlap fraction and the WHAM decorrelation time for the four US analogs, so that the choice of `k = 7.17 kcal mol⁻¹ Å⁻²` at `Δz = 1.0 Å` is either justified quantitatively or acknowledged as strong.
>
> **Inputs.**
> - US data: `../data/SUPP_umbrella_sampling/{resid}`
> - Analog subset: CYS⁰, SER, TYR⁰, ASP⁻ (four analogs).
> - Force constant and spacing used: `k = 7.17 kcal mol⁻¹ Å⁻²` (`3000 kJ mol⁻¹ nm⁻²`), `Δz = 1.0 Å`.
> - Code to produce the PMF : `../../supp_files/method_script/5-postprocessing/SUPP_umbrella_sampling/analysis/`.
>
> **Method.**
> 1. Identify the histogram data in the data directory.
> 2. **WHAM decorrelation.** For each window's `z(t)` in the production interval, compute the integrated autocorrelation time `τ_int` via the standard `[[ACF ESTIMATOR: e.g. emcee.autocorr.integrated_time or Sokal window]]`. Effective sample size per window: `N_eff = 150 ns / (2 τ_int)`.
> 3. Also compute the **theoretical thermal amplitude** for comparison: `σ_theo = √(k_B T / k) ≈ [[TO COMPUTE]] Å` at `T = 303.15 K` and quote the observed per-window standard deviation `σ_obs`. If `σ_obs > Δz`, adjacent windows are known to overlap well; if `σ_obs ≪ Δz`, the force constant is strong.
>
> **Outputs.**
> - `rev2/us_diagnostics/overlap_per_pair.csv` — columns `analog, window_i, window_i+1, BC`.
> - `rev2/us_diagnostics/per_analog_summary.csv` — one row per analog: `min_BC, mean_BC, max_BC, sigma_theo_A, mean_sigma_obs_A, mean_tau_int_ns, min_N_eff_per_window`.
> - `rev2/us_diagnostics/figures/BC_vs_window_<ANALOG>.pdf` — per-analog bar plot of overlaps.
> - `rev2/us_diagnostics/figures/tau_int_vs_window_<ANALOG>.pdf` — per-analog `τ_int` per window.
> - `rev2/us_diagnostics/REPORT.md` — one-paragraph summary phrased so it can be pasted into a new sub-block of Section S2.3 titled *Choice of US window spacing and force constant*.

---

## Analysis A5 — Reversibility: bulk↔bilayer transitions and residence times

**TODO2 link:** §1.5 (reversibility of "spontaneous partitioning").

> **Objective.** Quantify per-analog how many independent bulk-water↔bilayer transitions occur during production, so that the wording change to "quasi-equilibrium partitioning" for hydrophobic/aromatic analogs is backed by counts. The current text edit acknowledges qualitatively that (a) hydrophobic/aromatic analogs show few or no dissociations after insertion and (b) polar/charged analogs show many crossings; this analysis produces the numbers.
>
> **Inputs.** Same as A2.
>
> **Method.**
> 1. For each analog molecule and each replica, compute the per-frame folded depth `|z(t)|` on the full production trajectory (final 600 ns, or full 1 µs — state which is used; default: full 1 µs so the analysis captures early insertion events too).
> 2. Define **states** by the Region-I/II/III boundary vs Region-IV boundary: `bilayer = (|z| ≤ 30 Å)` (Regions I–III), `bulk = (|z| > 30 Å)` (Region IV). Use a hysteresis buffer of `[[HYSTERESIS: e.g. ±1 Å]]` on the 30 Å threshold to suppress spurious flicker crossings.
> 3. Count per-molecule **transitions**:
>    - `n_enter` = number of `bulk → bilayer` crossings,
>    - `n_leave` = number of `bilayer → bulk` crossings,
>    - `n_roundtrip` = min(`n_enter`, `n_leave`).
> 4. Aggregate per analog per replica by summing over the analog molecules present in that system (26 molecules at 0.2 M).
> 5. Compute **mean residence times** in each state per molecule: mean duration of a continuous `bilayer` interval; mean duration of a continuous `bulk` interval. Do not include the first and last intervals of each trajectory (their durations are censored).
> 6. Report per-analog **effective number of independent round-trips per trajectory** `N_ind = n_roundtrip / 26` averaged across replicas.
> 7. Mean/SE across the **three replicas** (not the 9 blocks — round-trips are trajectory-level quantities and 200-ns blocks are too short for many analogs). Use `s / √3`, `ddof = 1`.
> 8. Classify each analog:
>    - **reversibly sampled** if `N_ind ≥ 3` per replica,
>    - **quasi-equilibrium** if `0 < N_ind < 3`,
>    - **one-way insertion** if `n_leave = 0` in all three replicas.
>
> **Outputs.**
> - `rev2/reversibility/roundtrips_per_analog.csv` — one row per analog: `analog, n_enter_mean, n_enter_se, n_leave_mean, n_leave_se, N_ind_per_replica_mean, N_ind_per_replica_se, mean_residence_bilayer_ns, mean_residence_bulk_ns, class`.
> - `rev2/reversibility/figures/state_trace_<ANALOG>.pdf` — for a subset `[[SUBSET FOR TRACES: e.g. LEU, CYS0, SER, TYR0, ARG+, ASP-]]` show the state trace vs time for one analog molecule per replica.
> - `rev2/reversibility/REPORT.md` summarising the class of each analog and providing the paragraph text that will be dropped into the third Discussion limitation.

---

## Analysis A6 — Bilayer perturbation deltas at 0.2 M vs 0.1 M (compact table)

**TODO2 link:** §1.6.

> **Objective.** Assemble a compact main-text table showing the residual bilayer perturbation at 0.1 M vs 0.2 M for the eight residues that were resimulated at half concentration (VAL, ILE, LEU, MET, PHE, TYR⁰, TRP, ARG⁰), so Reviewer 2 can see that the perturbation is halved but not eliminated.
>
> **Inputs.**
> - Bilayer property outputs already used to make FigureS1–S6:
>   - phosphate-to-phosphate thickness `D_PP` (Å) per system
>   - area per lipid `A_PL` (Å²) per system
>   - deuterium order parameter profiles `S_CD` per system
>   - POPC-only reference values
> - Paths: `../data/SUPP_membrane_param/*csv`
>
> **Method.**
> 1. For each of the eight analogs and each of the two concentrations, compute the **paired absolute deviation** from the POPC-only reference:
>    - `|Δ D_PP| = D_PP^{system} − D_PP^{POPC}`
>    - `|Δ A_PL| = A_PL^{system} − A_PL^{POPC}`
>    - `%Δ S_CD = 100 · Σ_i |S_CD^{system}(i) − S_CD^{POPC}(i)| / Σ_i S_CD^{POPC}(i)` (matches the manuscript convention).
> 2. Mean and SE from the **paired 200-ns block-difference convention** (nine paired blocks) — do not use two-sample statistics.
> 3. Report a two-column table per property: `0.2 M mean ± SE` and `0.1 M mean ± SE`, plus a third column `ratio = 0.1M / 0.2M`.
>
> **Outputs.**
> - `rev2/concentration_deltas/perturbation_deltas.csv` — 8 rows (one per analog) × columns `dPP_02M_mean, dPP_02M_se, dPP_01M_mean, dPP_01M_se, dPP_ratio, dAPL_02M_mean, …, pctSCD_02M_mean, …, pctSCD_ratio`.
> - `rev2/concentration_deltas/TABLE_MAIN.tex` — a compact LaTeX table (`tabularx`, ≤ 5 columns wide) suitable for direct insertion into the main-text Results/Discussion of `Submission 3/manuscript/main.tex`.
> - `rev2/concentration_deltas/REPORT.md` — a short paragraph naming the property with the largest residual perturbation at 0.1 M and stating the ratio, ready to be dropped into the Discussion.

---

## Analysis A7 — Wimley–White octanol scale correlation panel

**TODO2 link:** §1.7 (octanol correlation).

> **Objective.** Add a scatter panel comparing the PMF-derived hydrocarbon-core transfer free energies to the Wimley–Creamer–White 1996 octanol side-chain transfer free energies, so Reviewer 2's point that the octanol scale is a canonical reference is addressed with data, not just wording.
>
> **Inputs.**
> - Computed PMFs (per analog, mean and SE on the 9-block convention): `../data/pmf_data/total/{resid}/`
> - Experimental octanol scale (Wimley, Creamer & White, *Biochemistry* 35, 5109–5124, 1996), digitised into a CSV: `../data/SUPP_hydrophobicity/oct_wat_data.dat`
>   - You must provide this file. Include the standard 20 amino acids + backbone contribution.
> - Analog-to-residue mapping: `Standard analogs`.
>
> **Method.**
> 1. For each analog with a matching residue in the octanol scale, define the computed transfer free energy as the Boltzmann-weighted regional free energy over the **hydrocarbon core** (`|z| ≤ 5 Å`), matching the convention already used in Section S16 of the supplemental. 
> Do the same Figure as in FigureS16_plot.ipynb but adding the octanol test as a third (pH=9) and fourth (pH=1) plot in the same figure.
> 4. **Do not** discard analogs from the fit for outlier reasons; label them in the plot.
>
> **Outputs.**
> - `rev2/octanol/octanol_vs_pmf.csv` with columns `residue, PMF_core_mean_kJmol, PMF_core_se_kJmol, exp_octanol_kJmol, exp_octanol_se_kJmol`.
> - `rev2/octanol/figures/FigureSx_octanol.png` — scatter with vertical + horizontal error bars, OLS line, `y = x` diagonal, and residue one-letter labels. Style consistent with the current `FigureS16.png`.
> - `rev2/octanol/REPORT.md` — one paragraph stating `R`, `ρ`, bootstrap CI, and the residues that deviate most, phrased for insertion into Section S6 of the supplemental.

---

## Analysis A8 — Full-population vs contact-free-subset PMF gap table (for the eight aggregating analogs)

**TODO2 link:** §1.2 (aggregation contribution to MacCallum disagreement).

> **Objective.** Quantify by how much the PMF changes when the associated (multi-analog contact) population is removed. This bounds the aggregation contribution to the MacCallum disagreement without adding a new figure to the supplemental.
>
> **Inputs.** Same as A2, plus per-frame inter-analog heavy-atom contact information (`d_min < 4.5 Å` criterion used in the manuscript).
>
> **Method.**
> 1. For each of the eight aggregating analogs (VAL, LEU, ILE, MET, PHE, TYR⁰, TRP, ARG⁰), partition every frame into `contact-free` and `associated` per analog molecule using the 4.5 Å criterion already implemented for Table 1 (Table S2/S3 in the current supplemental).
> 2. Compute the **contact-free-subset PMF** using only the contact-free frames as input to the density histogram, and the **full-population PMF** from all frames. Both use the same shared PMF convention.
> 3. Report the **residual PMF gap**:
>    - `Δ_min = PMF_full(z_min) − PMF_contact-free(z_min)` at the interfacial minimum depth,
>    - `RMS_gap = √⟨(PMF_full(z) − PMF_contact-free(z))²⟩` over the well-sampled depth range.
> 4. Mean and SE on the standard 9-block convention.
>
> **Outputs.**
> - `rev2/aggregation_gap/contact_free_gap.csv` — 8 rows: `analog, z_min_A, dGmin_full_kJmol, dGmin_cf_kJmol, delta_min_mean, delta_min_se, rms_gap_mean, rms_gap_se`.
> - `rev2/aggregation_gap/REPORT.md` — one-paragraph summary explicitly stating the mean absolute residual gap across the eight analogs and phrased for insertion into the Section S4 decomposition paragraph as evidence bounding the aggregation contribution.

## Reminders for whichever Claude runs the code

- **Do not invent input paths.** If a placeholder is left in the prompt, stop and ask.
- **Modify only one file** You're just authorized to modify the Figures_extra.ipynb notebook, each figure having independent cell.
- **Do not silently change the mean/SE convention.** If the analysis intrinsically cannot use the 9-block or 3-block convention, state the deviation in `REPORT.md` and use `ddof = 1` throughout.
- **Do not modify** any file under `Submission 3/`. Write all outputs into `rev2/<analysis>/`.
- **Do not add new figures or tables to the manuscript sources.** That is a separate step performed by a human (me) after review of the outputs.

DO NOT DO THE INSTRUCTIONS BELOW UNLESS YOU ARE IN sc-article-biophysj/ DIRECTORY. IF YOU ARE NOT SURE, ASK PERMISSION !!

## Post-analysis pass (do not do until A1–A8 outputs are in)

After the analyses above are complete, a further **text-only** pass will:
- Insert the actual numerical summaries from A1, A5, A6, A7, A8 into the Discussion paragraphs and the Section S4 decomposition already rewritten during the 2026-09-18 text-only pass.
- Add per-analog convergence flags (A1) and reversibility class labels (A5) to the PMF figure captions.
- Add the A6 table and A7 correlation panel to the manuscript / supplemental at that point.

Nothing above requires changing the current figure numbering. Any new figure introduced by A7 slots in as `FigureSx` (numbering to be finalised when the manuscript is switched to `elsarticle.cls` for BBA Biomembranes).
