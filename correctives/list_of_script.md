# List of scripts to correct

This file enumerates every script in the repository that still contains
hard-coded absolute paths, the `homoPOPC-<aa>-N` branch, the legacy
`traj=(4 ...)` / `traj=(4 5 6)` remappings, or a truncated/debug analog
list. Apply the canonical snippets below to every script in section
**"Scripts to fix"**.
build : add submit.sh (submit-local.sh and submit-alliancecan.sh)
---

## Canonical snippets to copy-paste

### Bash analog list (31 systems: 28 SC + 2 at 0.1 M + NONE = pure POPC)

```bash
aafile=( \
  "SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCN" "SCQ" \
  "SCF" "SCY" "SCW" "SCP" "GLYD" \
  "SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN" \
  "SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR" \
  "SCRN-1" "SCW-1" \
  "NONE" \
)

### Python analog list (31 systems: 28 SC + 2 at 0.1 M + NONE = pure POPC)

resname   = [ 
    "SCA", "SCV", "SCL", "SCI", "SCC", "SCM", "SCS", "SCT", "SCN", "SCQ",
    "SCF", "SCY", "SCW", "SCP", "GLYD",
    "SCHE", "SCHD", "SCDN", "SCEN", "SCKN", "SCRN",
    "SCHP", "SCD", "SCE", "SCCM", "SCYM", "SCK", "SCR",
    "SCRN-1", "SCW-1",
    "NONE"
    ] 

analogs = [
    "sca", "scv", "scl", "sci", "scc", "scm", "scs", "sct", "scn", "scq",
    "scf", "scy", "scw", "scp", "glyd",
    "sche", "schd", "scdn", "scen", "sckn", "scrn",
    "schp", "scd", "sce", "sccm", "scym", "sck", "scr",
    "scrn-1", "scw-1",
    "none",            # pure POPC reference
]
# Reviewer-facing layout. ${aa^^} is the upper-case label (e.g. SCV, SCRN, SCW-1),
# ${aa,,} is the lower-case label (e.g. scv, scrn, scw-1).
# All trajectories are hosted on Borealis (DOI: unavailable for now), one
# top-level directory per system. Each section file is compressed (.dcd.gz)
# and covers 100 ns (not 200 ns); there is no "ns" suffix in section names.
DIR="POPC-${aa^^}"

# PSF / DCD layout under $DIR :
#   $DIR/popc-${aa,,}.psf
#   $DIR/trajectory1/section401-500.dcd.gz
#   $DIR/trajectory1/section501-600.dcd.gz
#   $DIR/trajectory1/section601-700.dcd.gz
#   $DIR/trajectory1/section701-800.dcd.gz
#   $DIR/trajectory1/section801-900.dcd.gz
#   $DIR/trajectory1/section901-1000.dcd.gz
#   $DIR/trajectory2/...   (same 6 compressed sections)
#   $DIR/trajectory3/...   (same 6 compressed sections)
#
# In bash loops, decompress on the fly, e.g.:
#   for s in 401-500 501-600 601-700 701-800 801-900 901-1000; do
#     zcat "$DIR/trajectory${t}/section${s}.dcd.gz" > "section${s}.dcd"
#   done
# or feed `catdcd` directly via process substitution / a temp dir.

traj=(1 2 3)