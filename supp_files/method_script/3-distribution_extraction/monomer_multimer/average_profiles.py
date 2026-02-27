#!/usr/bin/env python3
import pandas as pd
import numpy as np

def load_and_normalize(filename):
    """
    Charge un fichier contenant deux colonnes :
        z   value
    Normalise value par sa somme totale.
    Retourne un DataFrame {z, norm}.
    """
    # Lecture : 2 colonnes seulement, séparateurs quelconques
    df = pd.read_csv(filename, sep=r"\s+", header=None, names=["z", "val"], usecols=[0, 1])

    # Normalisation (val/somme)
    total = df["val"].sum()
    # Si somme = 0 → laisser les valeurs intactes (cas rare)
    df["norm"] = df["val"] / total if total != 0 else df["val"]

    # On ne garde que z, norm
    df = df[["z", "norm"]]
    return df


def average_profiles(files):
    """Fusionne tous les DataFrames sur z, remplit vide par 0, moyenne."""
    merged = None

    for i, f in enumerate(files, start=1):
#        if i==1:
#            print('First file loaded :', f)
        df = load_and_normalize(f)
        df = df.rename(columns={"norm": f"norm_{i}"})

        if merged is None:
            merged = df
        else:
            # Fusion sur z, en gardant toutes les valeurs de z
            merged = pd.merge(merged, df, on="z", how="outer")
            
    # Remplir les z manquants par 0 (densité 0 dans ce fichier)
    merged = merged.fillna(0)

    # Moyenne des colonnes norm_i
    norm_cols = [c for c in merged.columns if c.startswith("norm_")]
    #print(norm_cols)
    merged["avg_norm"] = merged[norm_cols].mean(axis=1)
    #print(merged)
    # Tri en z
    merged = merged.sort_values("z").reset_index(drop=True)
    return merged


if __name__ == "__main__":
    import sys
    files = sys.argv[1:]
    if not files:
        print("Usage : python average_profiles.py fichier1.dat fichier2.dat ...")
        exit(1)

    result = average_profiles(files)
    print(result[["z", "avg_norm"]].to_string(index=False, header=False))

