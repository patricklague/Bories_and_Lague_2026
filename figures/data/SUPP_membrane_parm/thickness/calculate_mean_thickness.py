import pandas as pd
import numpy as np
import glob

# Liste des fichiers système
file_list = sorted(glob.glob("*thickness.dat"))

# Liste pour stocker les résultats
results = []

for filename in file_list:
    # Nom du système sans l'extension
    systeme_name = filename.split("-")[0]
    if filename.split("-")[1]=='1':
    	systeme_name = systeme_name+'-1'

    # Lecture du fichier
    df = pd.read_csv(filename, sep='\s+')

    # Filtrage des sections entre 400 et 1000
    #df_filtered = df[(df["#section"] >= 400) & (df["#section"] <= 1000)]
    #La commande au dessus a été mis en commentaire car les thickness ne contiennent que les sections 400 à 1000 par défaut.
    df_filtered=df.copy()
    # Récupérer toutes les colonnes "thickness-*"
    thickness_cols = [col for col in df_filtered.columns if col.startswith("thickness-")]

    # Calcul de la moyenne par batch de 20000 frames pour chaque trajectoire
    batches = [(0, 19999), (20000, 39999), (40000, 59999)]
    batch_means = []
    for start, end in batches:
        df_batch = df_filtered[(df_filtered["#frame"] >= start) & (df_filtered["#frame"] <= end)]
        for col in thickness_cols:
            vals = df_batch[col].dropna().values
            if len(vals) > 0:
                batch_means.append(np.mean(vals))

    batch_means = np.array(batch_means)
    mean_thickness = np.mean(batch_means)
    std_error = np.std(batch_means, ddof=1) / np.sqrt(len(batch_means))

    # Ajouter à la liste
    results.append([systeme_name, mean_thickness, std_error])

# Création du DataFrame final
summary_df = pd.DataFrame(results, columns=["name", "thickness", "std_error"])

# Sauvegarde
summary_df.to_csv("all_thicknesses.dat", sep="\t", index=False, float_format="%.4f")

