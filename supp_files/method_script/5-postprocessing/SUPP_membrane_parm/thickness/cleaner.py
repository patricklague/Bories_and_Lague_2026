import pandas as pd
import glob

# Liste des fichiers
file_list = sorted(glob.glob("scym-thickness.dat"))

for filename in file_list:
    df = pd.read_csv(filename, sep='\s+')

    # Retirer les frames 1 à 100
    df = df[df["#frame"] > 100].copy()

    # Réindexer les frames en soustrayant 100
    df["#frame"] = df["#frame"] - 101

    # Sauvegarder (écrase le fichier original)
    df.to_csv(filename, sep="\t", index=False, float_format="%.2f")