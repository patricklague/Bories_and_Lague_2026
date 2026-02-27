# Method Scripts

This directory contains all the scripts used to generate, analyze, and process the simulation data.
The workflow is divided into four sequential steps, each in its own subdirectory.

---

## 1-system_generation/

**Purpose:** Build the bilayer–solute systems (PSF/PDB files) for NAMD simulations.

**Entry point:** `build/launcher-run.sh`

`launcher-run.sh` loops over all amino acid analogs and calls `run.sh` for each one. `run.sh` orchestrates the full build in three steps:
1. Distribute solutes in the simulation box using `packmol-homoPOPC.sh` (with solute PDB files from `solutes/`)
2. Combine solutes with the bilayer system using `psfgen.inp` (VMD/psfgen)
3. Generate the extra-bonds file using `extraBondsFile.py`

**Contents:**
- `build/` — Build scripts (`launcher-run.sh`, `run.sh`, `packmol-homoPOPC.sh`, `psfgen.inp`, `extraBondsFile.py`, `update_step7.sh`)
- `solutes/` — PDB files for all amino acid side-chain analogs and associated topology/parameter files
- `charmm-gui.tgz` — Base POPC bilayer system from CHARMM-GUI

---

## 2-system_analysis/

**Purpose:** Analyze the MD trajectories: cell dimensions, membrane thickness, density profiles, and lipid order parameters.

**Entry point:** `analysis-run.sh`

`analysis-run.sh` loops over all analogs and trajectories, substituting placeholders in `analysis-self.sh` and running it. `analysis-self.sh` performs:
1. Trajectory centering (using `indexNoWat.vmd`, `center.vmd`, `indexNoMemb.vmd`)
2. Cell dimension extraction (plotted with `cell.gnu`)
3. Membrane thickness (using `thickness.vmd`, plotted with `thickness.gnu`)
4. POPC component density profiles (using `densityProfiles-popc.vmd` or `densityProfiles-water-total.vmd`)
5. Amino acid density profiles (using `densityProfiles-aa-total.vmd`)
6. Lipid acyl chain order parameters (using `orders-popc.vmd`)

The three density profile VMD scripts all depend on `density_profile.tcl` (VMD Density Profile Tool plugin, sourced at runtime).

---

## 3-distribution_extraction/

**Purpose:** Extract solute density distributions from the centered trajectories, separating monomers from multimers.

Contains two subdirectories:

### total/

**Entry point:** `analysis-run.sh`

Extracts the total (unseparated) amino acid density profiles. `analysis-run.sh` calls `analysis-self.sh` for each analog/trajectory, which uses `indexNoWat.vmd`, `center.vmd`, `indexNoMemb.vmd`, and `densityProfiles-aa-total.vmd` (+ `density_profile.tcl`).

### monomer_multimer/

**Entry point:** `launch_analysis_batched.sh`

Processes trajectories in batches to manage memory. For each batch, the script:
1. Removes water/ions (`indexNoWat.vmd`)
2. Centers the bilayer (`center.vmd`)
3. Removes the membrane (`indexNoMemb.vmd`)
4. Computes per-frame monomer vs. multimer density profiles (`densityProfiles-aa.vmd`, which sources `density_profile.tcl`)
5. Averages per-batch profiles (`average_profiles.py`)
6. Combines monomer and multimer profiles (`combine_profiles.py`)

---

## 4-pmf_calculation/

**Purpose:** Compute potentials of mean force (PMFs) from the density distributions extracted in step 3.

**Entry points:**
- `run.sh` — Computes PMFs from **total** distributions (calls `pmf-from-distribution-total.py`)
- `run-mono.sh` — Computes PMFs from **monomer** distributions (calls `pmf-from-distribution-mono.py`)
- `run-multi.sh` — Computes PMFs from **multimer** distributions (calls `pmf-from-distribution-multi.py`)

Each launcher loops over all analogs, copies the relevant distribution files from the analysis results, runs the corresponding Python PMF script, and stores the output PMF data.
