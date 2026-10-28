"""Phenotype statistics: one-way ANOVA, Tukey HSD and Benjamini-Hochberg.

Inputs : data/tissue_weights.csv, data/blood_glucose.csv
Output : outputs/phenotype_stats_full.csv  (12 endpoints x all three contrasts)

Every endpoint is reported with all three pairwise contrasts, including the
high-AGE versus low-AGE comparison, whether or not it is significant.
"""
import pandas as pd

from common import GROUPS, PAIRS, add_bh, anova_tukey, ensure_out_dir, load

ENDPOINTS = {
    "Body weight (g)": "body_weight_g",
    "SOL (mg)": "soleus_mg",
    "GAS (mg)": "gastrocnemius_mg",
    "GAS dry (mg)": "gastrocnemius_dry_mg",
    "EDL (mg)": "edl_mg",
    "Epididymal fat (mg)": "epididymal_fat_mg",
    "SOL/BW (mg/g)": "soleus_per_bw_mg_g",
    "GAS/BW (mg/g)": "gastrocnemius_per_bw_mg_g",
    "GAS dry/BW (mg/g)": "gastrocnemius_dry_per_bw_mg_g",
    "EDL/BW (mg/g)": "edl_per_bw_mg_g",
    "Fat/BW (mg/g)": "epididymal_fat_per_bw_mg_g",
}


def main():
    tissue = load("tissue_weights.csv")
    glucose = load("blood_glucose.csv")
    glucose = glucose.rename(columns={"blood_glucose_mg_dl": "value"})

    rows = []
    for label, column in ENDPOINTS.items():
        result = anova_tukey(tissue, column)
        result["metric"] = label
        rows.append(result)
    result = anova_tukey(glucose, "value")
    result["metric"] = "Blood glucose (mg/dL)"
    rows.append(result)

    frame = pd.DataFrame(rows).set_index("metric")
    frame = add_bh(frame)
    columns = (["n", "F", "ANOVA_p", "ANOVA_q_BH"]
               + [f"{k}_{a}_vs_{b}" for a, b in PAIRS
                  for k in ("Tukey_p", "q_BH", "diff", "ci_lo", "ci_hi")]
               + [f"{g}_{s}" for g in GROUPS for s in ("mean", "sd")])
    out = ensure_out_dir() / "phenotype_stats_full.csv"
    frame[columns].round(5).to_csv(out)
    print(f"wrote {out} ({frame.shape[0]} endpoints)")
    print("ANOVA q < 0.05:", frame.index[frame.ANOVA_q_BH < 0.05].tolist())
    for a, b in PAIRS:
        n_sig = int((frame[f"Tukey_p_{a}_vs_{b}"] < 0.05).sum())
        n_q = int((frame[f"q_BH_{a}_vs_{b}"] < 0.05).sum())
        print(f"  {a} vs {b}: Tukey p < 0.05 for {n_sig}, q < 0.05 for {n_q}")


if __name__ == "__main__":
    main()
