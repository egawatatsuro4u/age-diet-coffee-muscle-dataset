"""Pathway over-representation against an explicitly stated background.

Inputs : data/pathway_membership.csv, data/metabolome_long.csv,
         outputs/metabolome_stats_full.csv
Output : outputs/pathway_ORA_recomputed.csv

The background is the 54 metabolites quantified in this study. The 568-entry
target library is deliberately NOT used as the background, because most of its
entries were not detected and would inflate every enrichment.
"""
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from common import OUT_DIR, PAIRS, ensure_out_dir, load


def ora(selected, background, pathways, label):
    n_background = len(background)
    n_selected = len(selected)
    rows = []
    for pathway, members in pathways.items():
        in_background = [c for c in members if c in background]
        k = len([c for c in selected if c in in_background])
        p = hypergeom.sf(k - 1, n_background, len(in_background), n_selected) if k > 0 else 1.0
        rows.append({"set": label, "pathway": pathway, "background_N": n_background,
                     "pathway_K": len(in_background), "selected_n": n_selected,
                     "hits_k": k, "expected": n_selected * len(in_background) / n_background,
                     "p_hypergeom": p})
    frame = pd.DataFrame(rows)
    frame["q_BH"] = multipletests(frame.p_hypergeom, method="fdr_bh")[1]
    return frame.sort_values("p_hypergeom")


def main():
    membership = load("pathway_membership.csv")
    pathways = (membership[membership.in_measured_panel]
                .groupby("pathway").compound.apply(list).to_dict())
    background = sorted(load("metabolome_long.csv").compound.unique())
    stats_frame = pd.read_csv(OUT_DIR / "metabolome_stats_full.csv", index_col=0)

    selections = {
        "ANOVA_p<0.05": stats_frame.index[stats_frame.ANOVA_p < 0.05].tolist(),
        "AGEs_vs_Coffee_p<0.05": stats_frame.index[
            stats_frame["Tukey_p_AGEs_vs_AGEs+Coffee"] < 0.05].tolist(),
        "Control_vs_Coffee_p<0.05": stats_frame.index[
            stats_frame["Tukey_p_Control_vs_AGEs+Coffee"] < 0.05].tolist(),
    }
    frame = pd.concat([ora(sel, background, pathways, label)
                       for label, sel in selections.items()])
    out = ensure_out_dir() / "pathway_ORA_recomputed.csv"
    frame.round(4).to_csv(out, index=False)
    print(f"wrote {out} ({len(frame)} tests over {len(pathways)} pathways)")
    print(f"background N = {len(background)}; smallest q = {frame.q_BH.min():.3f}; "
          f"significant after BH: {int((frame.q_BH < 0.05).sum())}")


if __name__ == "__main__":
    main()
