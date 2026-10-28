"""Figure 3: Ward clustering of the z-scored metabolite matrix.

Input  : data/metabolome_long.csv
Outputs: outputs/Fig3_metabolome_heatmap.{png,pdf}

Rows (metabolites) are clustered; columns (pooled samples) are kept in fixed
group order so that the panel shows the metabolite blocks and not an inferred
sample grouping.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from common import COLOURS, LABELS, SAMPLES, ensure_out_dir, load, metabolome_wide

SAMPLE_GROUP = {s: g for g, members in
                zip(["Control", "AGEs", "AGEs+Coffee"],
                    [SAMPLES[:3], SAMPLES[3:6], SAMPLES[6:]])
                for s in members}


def main():
    wide = metabolome_wide(load("metabolome_long.csv"))
    wide.index = wide.index.get_level_values("compound")
    z_scored = wide.sub(wide.mean(axis=1), axis=0).div(wide.std(axis=1, ddof=1), axis=0)
    column_colours = [COLOURS[SAMPLE_GROUP[s]] for s in z_scored.columns]

    grid = sns.clustermap(
        z_scored, row_cluster=True, col_cluster=False, method="ward", metric="euclidean",
        cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, figsize=(7.2, 11.0),
        col_colors=column_colours, yticklabels=True,
        xticklabels=list(z_scored.columns),
        cbar_kws={"orientation": "horizontal"}, dendrogram_ratio=(0.11, 0.02),
        colors_ratio=0.018, cbar_pos=(0.30, 0.008, 0.18, 0.008))
    grid.ax_cbar.set_xlabel("z-score", fontsize=7, labelpad=2)
    grid.ax_cbar.tick_params(labelsize=6.5, length=2, pad=1.5)
    grid.ax_heatmap.set_xlabel("")
    grid.ax_heatmap.set_ylabel("")
    grid.ax_heatmap.tick_params(axis="y", labelsize=6.5, length=0)
    grid.ax_heatmap.tick_params(axis="x", labelsize=8, length=0)
    handles = [plt.Line2D([0], [0], marker="s", lw=0, markersize=8, color=c, label=g)
               for g, c in COLOURS.items()]
    grid.ax_col_dendrogram.legend(handles=handles, frameon=False, ncol=3, loc="center",
                                  fontsize=8)
    out_dir = ensure_out_dir()
    for suffix in ("png", "pdf"):
        grid.savefig(out_dir / f"Fig3_metabolome_heatmap.{suffix}", dpi=600,
                     bbox_inches="tight")
    plt.close(grid.figure)
    print(f"wrote Fig3 ({z_scored.shape[0]} metabolites x {z_scored.shape[1]} pooled samples)")


if __name__ == "__main__":
    main()
