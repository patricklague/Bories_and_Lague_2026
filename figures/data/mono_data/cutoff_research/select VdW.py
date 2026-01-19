from pymol import cmd
import math

# =====================
# Paramètres
# =====================
cutoff_vdw = 4.0
SC = "SCY"

cmd.select("vdw_pairs", "none")

# =====================
# Atomes lourds par SCY
# =====================
heavy_atoms = {}

for i in range(1, 27):
    sel = f"resname {SC} and resid {i} and not elem H"
    model = cmd.get_model(sel)
    if model.atom:
        heavy_atoms[i] = model.atom

# =====================
# Détection des paires VdW
# =====================
for i in heavy_atoms:
    for j in heavy_atoms:
        if j <= i:
            continue

        found = False
        for ai in heavy_atoms[i]:
            for aj in heavy_atoms[j]:
                dx = ai.coord[0] - aj.coord[0]
                dy = ai.coord[1] - aj.coord[1]
                dz = ai.coord[2] - aj.coord[2]
                if (dx*dx + dy*dy + dz*dz) ** 0.5 <= cutoff_vdw:
                    found = True
                    break
            if found:
                break

        if found:
            cmd.select(
                "vdw_pairs",
                f"vdw_pairs or (resname {SC} and resid {i}+{j})"
            )

