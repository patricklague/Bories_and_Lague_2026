# Bories and Lagüe, 2025

Supporting data and scripts for the article.

---

## figures/

Data, scripts, and outputs for all figures in the article.

- **data/** — Raw and processed data used for plotting (PMF data, density profiles, membrane parameters, aromatics orientation, hydrophobicity, monomer analysis)
- **scripts/** — Jupyter notebooks that generate the figures (PMF plots, density profiles, pKa, aromatics orientation, membrane parameters, monomer analysis, hydrophobicity)
- **plot/** — Output figure files (PNG)
- **draw/** — Hand-drawn figures (Pages files for article and abstract figures, ring orientation diagrams)

---

## supp_files/

Supplementary files: simulation inputs, outputs, parameters, and method scripts.

- **input_systems/** — Input PSF/PDB files for each amino acid analog system (one subdirectory per analog: sca, scc, ..., scym) and the pure POPC bilayer
- **output_systems/** — Output trajectory frames for each system (one subdirectory per analog), with helper scripts to extract frames (`get_last_frames.sh`, `get_trajs.sh`, `get_frame.tcl`)
- **parameter_files/** — CHARMM and NAMD topology/parameter files for both standard lipids (`toppar-charmm/`, `toppar-namd/`) and side-chain analogs (`sidechain-toppar-charmm/`, `sidechain-toppar-namd/`)
- **popc/** — Pure POPC reference data: density distributions and PMFs
- **method_script/** — All analysis scripts organized in four sequential steps (see [method_script/readme.md](supp_files/method_script/readme.md)):
  1. `1-system_generation/` — Build bilayer–solute systems
  2. `2-system_analysis/` — Trajectory analysis (centering, cell dimensions, thickness, density profiles, order parameters)
  3. `3-distribution_extraction/` — Extract solute density distributions (total and monomer/multimer)
  4. `4-pmf_calculation/` — Compute potentials of mean force from distributions