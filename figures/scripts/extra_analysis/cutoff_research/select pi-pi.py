from pymol import cmd
import math

# =====================
# Paramètres
# =====================
cutoff_pipi = 5.0
SC = "SCY"
aromatic_atoms = "name CG+CD1+CD2+CE1+CE2+CZ"

cmd.select("pipi", "none")

# =====================
# Calcul centres de masse aromatiques
# =====================
com = {}

for i in range(1, 27):
    sel = f"resname {SC} and resid {i} and ({aromatic_atoms})"
    if cmd.count_atoms(sel) == 0:
        continue
    com[i] = cmd.centerofmass(sel)

# =====================
# Détection des paires π–π
# =====================
for i in com:
    for j in com:
        if j <= i:
            continue

        dx = com[i][0] - com[j][0]
        dy = com[i][1] - com[j][1]
        dz = com[i][2] - com[j][2]
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)

        if dist <= cutoff_pipi:
            cmd.select("pipi", f"pipi or (resname {SC} and resid {i}+{j})")




