# Analysis code for the AGE-diet / coffee skeletal-muscle dataset

Code that reproduces every statistic and every figure reported in the Data in
Brief article *"Skeletal-muscle phenotype and GC/MS metabolome dataset from mice
fed a baked high-AGE diet with or without coffee intake."*

The data this code reads are the processed files of the companion deposit in the
Kyoto University Research Information Repository (KURENAI),
[TO BE COMPLETED: KURENAI handle and DOI of the deposit accompanying the revised
article]. That deposit supersedes an earlier one at
<http://hdl.handle.net/2433/301995> (doi:10.57723/kds627164), which remains
available under its own identifier. The raw GC/MS spectra are deposited in
MetaboLights under accession [TO BE COMPLETED: MetaboLights accession] and are
not required to run this code.

## Layout

```
analysis/
  common.py                 shared constants, loaders and the ANOVA/Tukey helper
  01_phenotype_stats.py     ANOVA, Tukey HSD and Benjamini-Hochberg, 12 endpoints
  02_metabolome_stats.py    the same per metabolite, plus the coffee-shift-index bootstrap
  03_pathway_ora.py         hypergeometric over-representation, explicit background
  04_figures.py             Figures 1, 2 and 4 and Supplementary Figure S1
  05_heatmap.py             Figure 3
data/                       processed CSV files, copied from the KURENAI deposit
outputs/                    everything the scripts write (not tracked in git)
```

## Reproducing the published values

```bash
conda env create -f environment.yml          # or: pip install -r requirements.txt
conda activate age-coffee-muscle
export PYTHONPATH=analysis
python analysis/01_phenotype_stats.py
python analysis/02_metabolome_stats.py
python analysis/03_pathway_ora.py            # needs the output of 02
python analysis/04_figures.py                # needs the outputs of 01-03
python analysis/05_heatmap.py
```

Scripts 01 to 03 write the four derived tables that are also in the deposit:
`phenotype_stats_full.csv`, `metabolome_stats_full.csv`,
`coffee_shift_index_bootstrap_CI.csv` and `pathway_ORA_recomputed.csv`. Running them
against the deposited input files reproduces the deposited tables exactly
(maximum absolute difference 0 across every numeric column).

The bootstrap is seeded (`BOOTSTRAP_SEED` in `common.py`), so the confidence
intervals are reproducible rather than merely similar. Changing the seed will
change the intervals in the last reported digits without changing any
conclusion.

## Things a reuser should know before running this

- **Metabolome samples are pools, not animals.** Each of the nine metabolome
  samples combines about 50 mg of left gastrocnemius from two animals, so that
  all six animals per group are represented. The sample identifiers
  (`L1`, `H1`, `HC1`, ...) are formally identical to the animal identifiers in
  `tissue_weights.csv` but denote different units. Joining the phenotype and
  metabolome files on those identifiers is wrong; use
  `metabolome_pool_map.csv`, which gives the pairing and the pool-level
  phenotype means.
- **No metabolite survives multiple-testing correction.** Five metabolites reach
  an uncorrected p < 0.05 and none reaches q < 0.05 (smallest q = 0.052). The
  metabolome analysis is hypothesis-generating.
- **The high-AGE and low-AGE groups do not differ** on any of the twelve
  phenotype endpoints or any of the 54 metabolites. Every significant difference
  in this dataset involves the coffee group.
- **The coffee shift index is unstable by construction.** Its denominator is the
  high-AGE minus low-AGE difference, which is not significant for any
  metabolite. The bootstrap interval has a median width of 22.4 and lies within
  [0, 1] for one metabolite out of 54. Treat the index as a direction label.
- **Group is confounded with injection order.** The nine samples were acquired
  in group order in a single batch, without randomisation and without pooled
  quality-control injections. `04_figures.py` writes the drift diagnostics.
- **Food and water intake are cage-level** with one cage per group, so they
  carry no within-group variance.

## Citing

Please cite the Data in Brief article and the KURENAI deposit. The archived
release of this repository has its own DOI; see the release notes.

## Licence

The code in this repository is released under the MIT Licence; see LICENSE.

The data files in `data/` and `outputs/` are part of the dataset deposit and
are released under CC BY 4.0, not under the MIT Licence. The two licences
differ because MIT is a software licence and CC BY 4.0 is the licence applied
to the deposited data; attribution is requested for both.
