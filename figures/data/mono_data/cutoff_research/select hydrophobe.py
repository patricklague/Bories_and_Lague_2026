from pymol import cmd
import math

# =====================
# Paramètres
# =====================
cutoff_hydrophobic = 4.5
SC = "SCY"

# Atomes carbones hydrophobes (exclut carbonyles, hétéroatomes)
hydro_sel = "elem C" #and not (name C=O+CO+O+N+S)"

cmd.select("hydrophobic_pairs", "none")

# =====================
# Récupération des atomes hydrophobes par SCY
# =====================
hydro_atoms = {}

for i in range(1, 27):
    sel = f"resname {SC} and resid {i} and ({hydro_sel})"
    model = cmd.get_model(sel)
    if len(model.atom) == 0:
        continue
    hydro_atoms[i] = model.atom

# =====================
# Détection des paires hydrophobes
# =====================
for i in hydro_atoms:
    for j in hydro_atoms:
        if j <= i:
            continue

        found = False
        for ai in hydro_atoms[i]:
            for aj in hydro_atoms[j]:
                dx = ai.coord[0] - aj.coord[0]
                dy = ai.coord[1] - aj.coord[1]
                dz = ai.coord[2] - aj.coord[2]
                dist = math.sqrt(dx*dx + dy*dy + dz*dz)

                if dist <= cutoff_hydrophobic:
                    found = True
                    break
            if found:
                break

        if found:
            cmd.select(
                "hydrophobic_pairs",
                f"hydrophobic_pairs or (resname {SC} and resid {i}+{j})"
            )


