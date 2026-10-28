"""Shared helpers for the analysis of the AGE-diet / coffee muscle dataset.

All scripts read from the deposited processed CSV files and write their outputs
to the directory given by OUT_DIR. No step depends on data that is not in the
deposit.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
OUT_DIR = Path(__file__).resolve().parents[1] / "outputs"

#: Group order used for every table and figure in the article.
GROUPS = ["Control", "AGEs", "AGEs+Coffee"]
#: Short labels used in the article text and figures.
LABELS = {"Control": "L", "AGEs": "H", "AGEs+Coffee": "H+C"}
#: Colours used consistently across all figures.
COLOURS = {"Control": "#4C72B0", "AGEs": "#C44E52", "AGEs+Coffee": "#55A868"}
#: The three pairwise contrasts, reported in full for every endpoint.
PAIRS = [("Control", "AGEs"), ("AGEs", "AGEs+Coffee"), ("Control", "AGEs+Coffee")]
#: Metabolome sample identifiers in injection order. Each is a pool of two
#: animals; see metabolome_pool_map.csv. They are NOT animal identifiers.
SAMPLES = ["L1", "L2", "L3", "H1", "H2", "H3", "HC1", "HC2", "HC3"]
#: Seed for the bootstrap, so that the deposited confidence intervals reproduce.
BOOTSTRAP_SEED = 20261014
BOOTSTRAP_REPLICATES = 10_000


def load(name):
    """Read a deposited CSV by file name."""
    return pd.read_csv(DATA_DIR / name)


def anova_tukey(df, value, group="group"):
    """One-way ANOVA plus Tukey HSD for all three pairwise contrasts.

    Returns a dict with the F statistic, the ANOVA p-value, the per-group n,
    mean and SD, and for each contrast the Tukey p-value, the mean difference
    and its 95% confidence interval.
    """
    arrays = [df.loc[df[group] == g, value].dropna().values for g in GROUPS]
    f_stat, p_value = stats.f_oneway(*arrays)
    tukey = pairwise_tukeyhsd(df[value].values, df[group].values)
    table = pd.DataFrame(tukey.summary().data[1:], columns=tukey.summary().data[0])
    out = {"F": f_stat, "ANOVA_p": p_value, "n": [len(a) for a in arrays]}
    for a, b in PAIRS:
        row = table[((table.group1 == a) & (table.group2 == b))
                    | ((table.group1 == b) & (table.group2 == a))].iloc[0]
        out[f"Tukey_p_{a}_vs_{b}"] = float(row["p-adj"])
        out[f"diff_{a}_vs_{b}"] = float(row["meandiff"])
        out[f"ci_lo_{a}_vs_{b}"] = float(row["lower"])
        out[f"ci_hi_{a}_vs_{b}"] = float(row["upper"])
    for g, arr in zip(GROUPS, arrays):
        out[f"{g}_mean"] = arr.mean()
        out[f"{g}_sd"] = arr.std(ddof=1)
    return out


def add_bh(frame):
    """Add Benjamini-Hochberg q-values across rows for ANOVA and each contrast.

    Correction is applied across the rows of `frame`, that is across endpoints
    or across metabolites, separately for the ANOVA and for each contrast.
    """
    frame["ANOVA_q_BH"] = multipletests(frame["ANOVA_p"], method="fdr_bh")[1]
    for a, b in PAIRS:
        frame[f"q_BH_{a}_vs_{b}"] = multipletests(frame[f"Tukey_p_{a}_vs_{b}"],
                                                  method="fdr_bh")[1]
    return frame


def metabolome_wide(long_frame):
    """Pivot the tidy metabolome table to metabolites x pooled samples."""
    wide = long_frame.pivot_table(index=["compound", "category"], columns="sample_id",
                                  values="peak_area_ratio_corrected")
    return wide[SAMPLES]


def ensure_out_dir():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUT_DIR
