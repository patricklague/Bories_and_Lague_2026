from pymol import cmd
import math

# =====================
# Paramètres
# =====================
cutoff_com = 7.0
SC = "SCY"

cmd.select("com_pairs", "none")

# =====================
# Calcul des centres de masse
# =====================
com = {}

for i in range(1, 27):
    sel = f"resname {SC} and resid {i}"
    if cmd.count_atoms(sel) == 0:
        continue
    com[i] = cmd.centerofmass(sel)

# =====================
# Détection des paires
# =====================
for i in com:
    for j in com:
        if j <= i:
            continue

        dx = com[i][0] - com[j][0]
        dy = com[i][1] - com[j][1]
        dz = com[i][2] - com[j][2]
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)

        if dist <= cutoff_com:
            cmd.select(
                "com_pairs",
                f"com_pairs or (resname {SC} and resid {i}+{j})"
            )

