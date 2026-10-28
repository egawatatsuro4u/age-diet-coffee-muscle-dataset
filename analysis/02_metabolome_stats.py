"""Metabolome statistics and the bootstrap interval for the coffee shift index.

Inputs : data/metabolome_long.csv
Outputs: outputs/metabolome_stats_full.csv      (54 metabolites x all contrasts)
         outputs/coffee_shift_index_bootstrap_CI.csv  (index, CI and instability)

The coffee shift index (H - HC)/(H - L) has the high-AGE minus low-AGE difference in
its denominator. That difference is not significant for any metabolite in this
dataset, so the index is reported only with its bootstrap interval and with
diagnostics of how often the denominator changes sign.
"""
import numpy as np
import pandas as pd

from common import (BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED, GROUPS, PAIRS, SAMPLES,
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


def coffee_shift_index_bootstrap(wide):
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    n_rep = BOOTSTRAP_REPLICATES
    rows = []
    for (compound, category), row in wide.iterrows():
        control, high, coffee = (row[SAMPLE_SETS[g]].values.astype(float) for g in GROUPS)
        index = (high.mean() - coffee.mean()) / (high.mean() - control.mean())
        draws = [arr[rng.integers(0, 3, (n_rep, 3))].mean(axis=1)
                 for arr in (control, high, coffee)]
        boot_control, boot_high, boot_coffee = draws
        denominator = boot_high - boot_control
        boot_index = (boot_high - boot_coffee) / denominator
        lo, hi = np.nanpercentile(boot_index, [2.5, 97.5])
        rows.append({
            "compound": compound,
            "category": category,
            "CoffeeShiftIndex": index,
            "boot_CI_lo": lo,
            "boot_CI_hi": hi,
            "CI_width": hi - lo,
            "denom_sign_flip_frac": float(np.mean(
                np.sign(denominator) != np.sign(high.mean() - control.mean()))),
            "frac_boot_outside_0_1": float(np.mean((boot_index < 0) | (boot_index > 1))),
        })
    return pd.DataFrame(rows).set_index("compound")


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
    shift_index = coffee_shift_index_bootstrap(wide)
    shift_index = shift_index.join(stats_frame[["ANOVA_p", "ANOVA_q_BH", "Tukey_p_Control_vs_AGEs"]])
    shift_index.round(4).to_csv(out_dir / "coffee_shift_index_bootstrap_CI.csv")

    n_anova = int((stats_frame.ANOVA_p < 0.05).sum())
    n_q = int((stats_frame.ANOVA_q_BH < 0.05).sum())
    print(f"metabolites with ANOVA p < 0.05: {n_anova}; surviving BH: {n_q} "
          f"(smallest q = {stats_frame.ANOVA_q_BH.min():.3f})")
    print(f"coffee shift index: median CI width {shift_index.CI_width.median():.1f}, "
          f"CI within [0, 1] for {int(((shift_index.boot_CI_lo > 0) & (shift_index.boot_CI_hi < 1)).sum())} "
          f"of {len(shift_index)} metabolites")


if __name__ == "__main__":
    main()
