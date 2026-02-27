from pymol import cmd

cutoff = 4.5
SC='SCY'
cmd.select("close", "none")

for i in range(1, 27):
    sc = f"sc{i}"
    cmd.select(sc, f"resname {SC} and resid {i} and not elem H")
    cmd.select(
        "close",
        f"close or (byres ((resname {SC} and not elem H) within {cutoff} of {sc}) and not {sc})"
    )

