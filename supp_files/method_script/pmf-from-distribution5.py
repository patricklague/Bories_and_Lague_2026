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
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=["x", "y"])
    all_x_values.update(abs(x) for x in df["x"])

x_full = sorted(all_x_values)
x_df_template = pd.DataFrame({"x": x_full})

# Liste des résultats à fusionner
all_results = [x_df_template.copy()]
curve_index = 1  # Compteur pour y1, y2, ..., y6

for filename in file_list:
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=["x", "y"])

    for sign, subset in [("neg", df[df["x"] < 0]), ("pos", df[df["x"] > 0])]:
        subset = subset.copy()
        subset["x"] = subset["x"].abs()

        label = f"y{curve_index}"
        raw_col = f"{label}_raw"
        norm_col = f"{label}_norm"
        pmf_col = f"PMF_{label}"

        # Grouper par x (si doublons), moyenne
        grouped = subset.groupby("x")["y"].mean().reset_index().rename(columns={"y": raw_col})
        grouped = pd.merge(x_df_template, grouped, on="x", how="left").fillna(0)

        # Normalisation
        total = grouped[raw_col].sum()
        grouped[norm_col] = grouped[raw_col] / total if total != 0 else grouped[raw_col]

        # PMF brut
        pmf_raw = -k_B * T * np.log(grouped[norm_col] + epsilon)

        # Recalage <PMF> entre x = 40 et x = 50
        mask_window = (grouped["x"] >= 40) & (grouped["x"] <= 50)
        pmf_avg = pmf_raw[mask_window].mean()
        grouped[pmf_col] = pmf_raw - pmf_avg

        # Conserver les colonnes
        final_df = grouped[["x", raw_col, norm_col, pmf_col]]
        all_results.append(final_df.drop(columns="x"))

        curve_index += 1

# Fusion finale
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
