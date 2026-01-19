from pymol import cmd
import math

# =====================
# Paramètres
# =====================
cutoff_hbond = 3.5
SC = "SCY"

cmd.select("hbond_pairs", "none")

# =====================
# Donneurs / accepteurs par SCY
# =====================
donors = {}
acceptors = {}

for i in range(1, 27):
    sel_d = f"resname {SC} and resid {i} and elem N+O"
    sel_a = f"resname {SC} and resid {i} and elem O+N"

    d_model = cmd.get_model(sel_d)
    a_model = cmd.get_model(sel_a)

    if d_model.atom and a_model.atom:
        donors[i] = d_model.atom
        acceptors[i] = a_model.atom

# =====================
# Détection des paires H-bond
# =====================
for i in donors:
    for j in acceptors:
        if j <= i:
            continue

        found = False

        # i donneur -> j accepteur
        for d in donors[i]:
            for a in acceptors[j]:
                dx = d.coord[0] - a.coord[0]
                dy = d.coord[1] - a.coord[1]
                dz = d.coord[2] - a.coord[2]
                if (dx*dx + dy*dy + dz*dz) ** 0.5 <= cutoff_hbond:
                    found = True
                    break
            if found:
                break

        # j donneur -> i accepteur
        if not found:
            for d in donors[j]:
                for a in acceptors[i]:
                    dx = d.coord[0] - a.coord[0]
                    dy = d.coord[1] - a.coord[1]
                    dz = d.coord[2] - a.coord[2]
                    if (dx*dx + dy*dy + dz*dz) ** 0.5 <= cutoff_hbond:
                        found = True
                        break
                if found:
                    break

        if found:
            cmd.select(
                "hbond_pairs",
                f"hbond_pairs or (resname {SC} and resid {i}+{j})"
            )
