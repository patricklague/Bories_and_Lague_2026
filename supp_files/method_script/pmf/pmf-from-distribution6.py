import pandas as pd
import numpy as np
import glob

# Constantes
k_B = 0.008314  # kJ/mol·K
T = 303.15
epsilon = 1e-5

# Lecture des fichiers
file_list = sorted(glob.glob("trajectory*.dat"))

# Construire la grille x positive commune
all_x_values = set()
for filename in file_list:
    df = pd.read_csv(filename, sep='\s+', header=None, names=["x", "y"])
    all_x_values.update(abs(x) for x in df["x"])

x_full = sorted(all_x_values)
x_df_template = pd.DataFrame({"x": x_full})

# Liste des résultats à fusionner
all_results = [x_df_template.copy()]
curve_index = 1  # Compteur pour y1, y2, ..., y6

for filename in file_list:
    df = pd.read_csv(filename, sep='\s+', header=None, names=["x", "y"])

    label = f"y{curve_index}"
    raw_col = f"{label}_raw"
    norm_col = f"{label}_norm"
    pmf_col = f"PMF_{label}"

    # --- 1. Regrouper tous les x (negatifs + positifs) et moyenner les doublons ---
    df_all = df.copy()
    df_all_grouped = df_all.groupby("x")["y"].mean().reset_index().rename(columns={"y": raw_col})

    # Aligner avec le template des x
    df_all_grouped = pd.merge(x_df_template, df_all_grouped, on="x", how="left").fillna(0)

    # --- 2. Normalisation sur tous les x ---
    total_all = df_all_grouped[raw_col].sum()
    df_all_grouped[norm_col] = df_all_grouped[raw_col] / total_all if total_all != 0 else df_all_grouped[raw_col]

    # --- 3. Séparer ensuite en x+ et x-, mais sur la distribution normalisée ---
    df_neg = df_all_grouped[df_all_grouped["x"] < 0].copy()
    df_pos = df_all_grouped[df_all_grouped["x"] > 0].copy()

    # Convertir les x en valeurs absolues
    df_neg["x"] = df_neg["x"].abs()
    df_pos["x"] = df_pos["x"].abs()

    # --- 4. Regrouper les deux moitiés en fonction de |x| et faire la moyenne ---
    # pour éviter des problèmes d’alignement, réindexer sur la grille x_df_template
    df_neg = pd.merge(x_df_template, df_neg[["x", norm_col]], on="x", how="left").fillna(0)
    df_pos = pd.merge(x_df_template, df_pos[["x", norm_col]], on="x", how="left").fillna(0)

    # moyenne entre distributions normalisées
    df_avg = x_df_template.copy()
    df_avg[norm_col] = (df_neg[norm_col] + df_pos[norm_col]) / 2

    # --- 5. Calcul de la PMF à partir de la distribution normalisée moyenne ---
    pmf_raw = -k_B * T * np.log(df_avg[norm_col] + epsilon)

    # Recalage (*) : <PMF> sur fenêtre 40–50
    mask_window = (df_avg["x"] >= 40) & (df_avg["x"] <= 50)
    pmf_avg = pmf_raw[mask_window].mean()

    df_avg[pmf_col] = pmf_raw - pmf_avg

    # --- 6. Stocker le résultat final ---
    final_df = df_avg[["x", norm_col, pmf_col]].copy()
    all_results.append(final_df.drop(columns="x"))

    curve_index += 1

# Fusion de toutes les colonnes finales
result = pd.concat(all_results, axis=1)

# Sauvegarde
result.to_csv("pmf.dat", sep="\t", index=False, float_format="%.6f")

# Moyennes et erreurs sur PMF
pmf_cols = [col for col in result.columns if col.startswith("PMF_y")]
pmf_array = result[pmf_cols].values

pmf_stats = result[["x"]].copy()
pmf_stats["PMF_mean"] = np.mean(pmf_array, axis=1)
pmf_stats["std_dev"] = np.std(pmf_array, axis=1, ddof=1)
pmf_stats["std_error"] = pmf_stats["std_dev"] / np.sqrt(len(pmf_cols))

pmf_stats.to_csv("pmf_moyen.dat", sep="\t", index=False, float_format="%.6f")

print("\nPMF moyen et erreurs :")
print(pmf_stats)
