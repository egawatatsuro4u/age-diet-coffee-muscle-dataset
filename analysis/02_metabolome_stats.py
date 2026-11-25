"""Metabolome statistics: ANOVA, Tukey HSD and Benjamini-Hochberg correction.

Inputs : data/metabolome_long.csv
Output : outputs/metabolome_stats_full.csv  (54 metabolites x all contrasts)

Treatment was assigned to one cage per group, so these model outputs are
conditional calculations that describe dispersion within a cage; they do not
establish a treatment effect.
"""
import numpy as np
import pandas as pd

from common import (GROUPS, PAIRS, SAMPLES,
                    add_bh, anova_tukey, ensure_out_dir, load, metabolome_wide)

SAMPLE_SETS = {"Control": SAMPLES[:3], "AGEs": SAMPLES[3:6], "AGEs+Coffee": SAMPLES[6:]}


def per_metabolite_stats(long_frame):
    rows = []
    for (compound, category), group in long_frame.groupby(["compound", "category"], sort=False):
        result = anova_tukey(group, "peak_area_ratio_corrected")
        result["compound"] = compound
        result["category"] = category
        rows.append(result)
    frame = pd.DataFrame(rows).set_index("compound")
    return add_bh(frame)


def main():
    long_frame = load("metabolome_long.csv")
    stats_frame = per_metabolite_stats(long_frame)
    columns = (["category", "n", "F", "ANOVA_p", "ANOVA_q_BH"]
               + [f"{k}_{a}_vs_{b}" for a, b in PAIRS
                  for k in ("Tukey_p", "q_BH", "diff", "ci_lo", "ci_hi")]
               + [f"{g}_{s}" for g in GROUPS for s in ("mean", "sd")])
    out_dir = ensure_out_dir()
    stats_frame[columns].round(6).to_csv(out_dir / "metabolome_stats_full.csv")

    wide = metabolome_wide(long_frame)

    n_anova = int((stats_frame.ANOVA_p < 0.05).sum())
    n_q = int((stats_frame.ANOVA_q_BH < 0.05).sum())
    print(f"metabolites with ANOVA p < 0.05: {n_anova}; surviving BH: {n_q} "
          f"(smallest q = {stats_frame.ANOVA_q_BH.min():.3f})")


if __name__ == "__main__":
    main()
