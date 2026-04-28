# Method Scripts

This directory contains all the scripts used to generate, analyze, and process the simulation data.
The workflow is divided into five sequential steps, each in its own subdirectory.

---

## 1-system_generation/

**Purpose:** Build the bilayer–solute systems (PSF/PDB files) for NAMD simulations.

**Entry point:** `build/run.sh`

`run.sh` loops over all amino acid analogs and calls `main.sh` for each one. `main.sh` orchestrates the full build in three steps:
1. Distribute solutes in the simulation box using `packmol-POPC.sh` (with solute PDB files from `solutes/`)
2. Combine solutes with the bilayer system using `psfgen.inp` (VMD/psfgen)
3. Generate the extra-bonds file using `extraBondsFile.py`

**Contents:**
- `build/` — Build scripts (`run.sh`, `main.sh`, `packmol-POPC.sh`, `psfgen.inp`, `extraBondsFile.py`, `update_step7.sh`)
- `solutes/` — PDB files for all amino acid side-chain analogs and associated topology/parameter files
- `charmm-gui.tgz` — Base POPC bilayer system from CHARMM-GUI

---

## 2-system_analysis/

**Purpose:** Analyze the MD trajectories: cell dimensions, membrane thickness, density profiles, and lipid order parameters.

**Entry point:** `analysis-run.sh`

`analysis-run.sh` loops over all analogs and trajectories, substituting placeholders in `analysis_per_system.sh` and running it. `analysis_per_system.sh` performs:
1. Trajectory centering (using `indexNoWat.vmd`, `center.vmd`)
2. Cell dimension extraction
3. Membrane thickness (using `thickness.vmd`)
4. POPC component density profiles (using `densityProfiles-popc.vmd` or `densityProfiles-water-total.vmd`)
5. Lipid acyl chain order parameters (using `orders-popc.vmd`)

Orientation for aromatics raw data extraction are contained in `aromatic_orientation/`

The density profile VMD scripts depend on `density_profile.tcl` (VMD Density Profile Tool plugin, sourced at runtime).


---

## 3-distribution_extraction/

**Purpose:** Extract per-frame solute positions from the centered trajectories and bin them into z-density distributions, identifying monomers from multimers via a contact cutoff.

**Entry point:** `counting_script.ipynb`

`counting_script.ipynb` is the master notebook that drives the whole step. It launches `get_trajectory.sh` for each analog/trajectory, which in turn:
1. Removes water/ions (`indexNoWat.vmd`)
2. Centers the bilayer (`center.vmd`)
3. Removes the membrane (`indexNoMemb.vmd`)
4. Writes the per-frame solute positions and per-frame contact counts at the 4.5 Å cutoff.

Once the raw contact files are produced, `raw_to_distribution.py` can be used to bin them into the per-trajectory z-density distributions consumed by step 4.

---

## 4-pmf_calculation/

**Purpose:** Compute potentials of mean force (PMFs) from the density distributions extracted in step 3.

**Entry point:** `run.sh`

`run.sh` loops over all analogs, copies the relevant distribution files from the analysis results, runs `pmf-from-distribution.py`, and stores the output PMF data. The script handles the `total`, `monomer_4.5A`, and `multimer_4.5A` distribution modes through its parameters.

---

## 5-data_formatting/

**Purpose:** Take the raw outputs of steps 1–4 and turn them into the artifacts that are actually shipped (Borealis trajectory bundle + the curated `figures/data/` tree). This step is the only one allowed to write outside `supp_files/method_script/`.

**Contents:**

- `borealis_formatting/` — Concatenate and compress the per-100 ns trajectory sections into the public Borealis layout (`POPC-<AA>/trajectory{1,2,3}/section{401-500,501-600,601-700,701-800,801-900,901-1000}.dcd.gz`) and copy the matching `popc-<aa>.psf` next to it.

- `SUPP_membrane_param_formatting/` — Per-system extraction and aggregation scripts that turn the raw `analysis_per_system.sh` outputs into the canonical `figures/data/SUPP_membrane_parm/` tables. Includes the `get_*.py` extractors (thickness, area per lipid, density profiles, deuterium order parameters), their `run.sh` drivers, and the `compute_*.py` aggregators that produce the `computed_*.csv` summaries.

- `aromatics_orientation_formatting/` — Convert the per-frame aromatic-ring atom coordinates produced in step 2 into the `θ` / `φ` ring-orientation angles consumed by `figures/scripts/aromatics_orientation_plots.ipynb`. Includes `Ring_xyz_preparation.ipynb`, `ring_orientation_analysis.py`, `freq_angle_analysis.py`, `extract_top_angles.sh`, and the `run.sh` driver.

All scripts in this step apply the canonical naming conventions defined in `correctives/list_of_script.md` and emit headers compatible with `correctives/fix_header.py` so that no post-hoc header rewriting is needed.
