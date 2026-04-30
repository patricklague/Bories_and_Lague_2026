import pandas as pd
import numpy as np
import glob

# Liste des fichiers système
file_list = sorted(glob.glob("../../../../../figures/data/SUPP_membrane_parm/area_per_lipid/*-apl.dat"))

# Liste pour stocker les résultats
results = []

for filename in file_list:
    # Nom du système sans l'extension
    systeme_name = filename.split("/")[-1].split("-")[0]
    if filename.split("-")[1] == "1":
    	systeme_name = systeme_name + "-1"

    # Lecture du fichier
    df = pd.read_csv(filename, sep='\s+')

    # Récupérer toutes les colonnes "apl*"
    apl_cols = [col for col in df.columns if col.startswith("apl")]

    # Calcul de la moyenne par batch de 200 sections pour chaque trajectoire
    batches = [(401, 600), (601, 800), (801, 1000)]
    batch_means = []
    for start, end in batches:
        df_batch = df[(df["#section"] >= start) & (df["#section"] <= end)]
        for col in apl_cols:
            vals = df_batch[col].dropna().values
            if len(vals) > 0:
                batch_means.append(np.mean(vals))

    batch_means = np.array(batch_means)
    mean_apl = np.mean(batch_means)
    std_error = np.std(batch_means, ddof=1) / np.sqrt(len(batch_means))

    # Ajouter à la liste
    results.append([systeme_name, mean_apl, std_error])

# Création du DataFrame final
summary_df = pd.DataFrame(results, columns=["name", "apl", "std_error"])

# Sauvegarde
summary_df.to_csv("../../../../../figures/data/SUPP_membrane_parm/area_per_lipid/all_apl.dat", sep="\t", index=False, float_format="%.4f")

