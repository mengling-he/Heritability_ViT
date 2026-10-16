# Heritability_ViT

Copy of `Heritability_Unet` with scripts switched to the **ViT** segmentation
features. Every script starts with

```r
segmentation_method <- "vit"  # Change to "unet" when needed
base_dir <- switch(segmentation_method,
                   "unet" = "Heritability_Unet",
                   "vit"  = "Heritability_ViT")
```

and all input/output file names are built with `sprintf("..._%s...", segmentation_method)`,
so setting `"unet"` reads/writes the `Heritability_Unet` folder instead.
ViT raw files: `raw_data/*_features_vit.csv` (date column `Analyse.Date`).

---

Original description (U-Net):

Heritability analysis of longitudinal poplar phenotypes derived from **U-Net
segmentation**. This is a simplified re-run of the `Heritability_update`
(ViT / thresholding) pipeline on the new U-Net feature data, dropping the DMD
and bootstrap analyses.

## Layout

```
raw_data/                       U-Net side & top feature CSVs + reference files
  *.side_features.csv, *.top_features.csv
  slow_grower.txt               reference slow-grower list (from ViT pipeline)
  h2_timefixed_day_threshold.csv reference thresholding-method H2 (for Fig 3)
Code/                           shared functions (heritability_fun.R, plot_fun.R, FPCA_code/)
Analysis/
  1_preprocessing/              combine side+top, slow growers, daily & date-group imputation
  2_heritability/               per-date and date-as-covariate H2; trajectory plot
  3_threshold_comparison/       U-Net vs thresholding-method and vs ViT H2
  4_FPCA/                       functional PCA scores + their H2 (traits z-scored first)
  5_LDA/                        LDA scores + H2 (incl. leave-one-genotype-out CV)
```
Each step folder has `script/` and `result/` (step 1 also has `data/`).

## Run order

Run from the **renv project root** (the `Heritability/` folder, which holds
`.Rprofile` + `renv.lock`) so the renv library activates:

```bash
cd /path/to/Heritability      # renv root
Rscript Heritability_ViT/Analysis/1_preprocessing/script/1_preprocess.R
Rscript Heritability_ViT/Analysis/2_heritability/script/2_heritability.R
Rscript Heritability_ViT/Analysis/3_threshold_comparison/script/3_compare_threshold.R
Rscript Heritability_ViT/Analysis/3_threshold_comparison/script/3c_recompute_vit_and_compare.R
Rscript Heritability_ViT/Analysis/4_FPCA/script/4_fpca.R
Rscript Heritability_ViT/Analysis/5_LDA/script/5_lda.R
```
Steps 2–5 depend on step 1; step 3 / 3c also read step 2; steps 4–5 read step 2.
Scripts locate their own paths, so the working directory only matters for renv.

## Models (Code/heritability_fun.R)

- per-date:           `trait ~ (1|Genotype)`                                   (`calc_H2_timeseries`)
- date-as-covariate:  `trait ~ Day + (1|Genotype) + (1|Genotype:Replicate)`    (`calc_H2_fixedtime`)
- FPC / single score: `score ~ (1|Genotype)`                                   (`calc_H2_1dim`)

`Genotype`/`Replicate` are parsed from `Plant.Info` (`^\d+_(.+)_(\d+)$`).

## Key outputs / report-figure mapping

| Report figure | File |
|---|---|
| Fig 2 missing pattern | `1_preprocessing/result/missing_pattern.png` |
| Fig 4 H2 trajectories (Side/Top) | `2_heritability/result/h2_series_day_vit_by_view.png` |
| H2 side-vs-top sorted by mean | `2_heritability/result/h2_timefixed_day_vit_by_view.png` |
| Fig 3 Threshold vs U-Net | `3_threshold_comparison/result/h2_compare_threshold_vit.png` |
| Fig 5 FPCA vs LMM | `4_FPCA/result/h2_FPCA_vs_LMM_vit.png` |
| Fig 6 LD1/LD2 loadings | `5_LDA/result/lda_loadings_LD{1,2}_vit.png` |

## Notes on this U-Net run

- 446 plants, 60 genotypes. **12** plants flagged as slow growers (side-area
  growth < 5%); the ViT reference list had 13 — the only difference is
  `63_GW-9762_3` (see `1_preprocessing/result/slow_grower_comparison.csv`).
  Viable plants used downstream: **434**.
- The daily grid uses the **31 observed dates** (the calendar days on which at
  least one plant was imaged), not a full day-by-day calendar. Each plant is
  filled in on every observed date within its measured span.
- Within-day aggregation uses the **median** ("middle point") of multiple images
  on the same day (max 2 images/day here, so median = midpoint).
- Daily alignment uses internal linear interpolation **plus end extrapolation**
  (`zoo::na.approx(..., rule = 2)`): leading/trailing gaps are filled by carrying
  the plant's first/last observed value. The imputed daily grid therefore has
  **no NA**, so all 31 observed dates (2024-05-26 to 2024-06-30) are complete.
- **Complete-date filtering:** the per-timepoint and LMM-with-time analyses keep
  only dates where *every* sample has a value. With end-extrapolation this now
  retains **all 31 dates**. Applied in step 2, the ViT recompute (step 3c), and
  the LDA H2 (step 5). **FPCA (step 4) is unaffected** — it uses the un-imputed
  combined trajectories. (Without extrapolation the filter had left only 16 dates,
  2024-06-03 to 2024-06-18, because one plant — 297_BESC-284_4 — was first imaged
  on 06-03 and two were last imaged on 06-18.)
- Note: step 3's comparison against the *precomputed* threshold / ViT references
  is this-pipeline-vs-old-pipeline; the apples-to-apples comparison is **step 3c**,
  where ViT is re-run through this pipeline (same median + extrapolation).
- Top raw-trait heritabilities (date-as-covariate): `top_saturation_mean` 0.48,
  `side_solidity` 0.46, `top_blue_yellow_mean` 0.43.
- Date-as-covariate H2 comparisons (step 3): vs **Threshold** roughly even
  (U-Net higher in 28/73 traits, mean diff +0.003); vs **ViT** the ViT pipeline
  is generally higher (U-Net higher in only 21/81 traits, mean diff −0.031).
  See `h2_compare_vit_unet.png` and the three-way `h2_compare_threshold_vit_unet.png`
  (ViT reference: `../Heritability_update/.../h2_timefixed_day_vit.csv`).
- FPC1 heritability exceeds the corresponding raw-trait H2 for the strongest
  traits (e.g. `top_saturation_mean` 0.68 vs 0.48).
- **FPCA inputs are z-scored** (step 4). Large-magnitude traits like `side_area`
  (~1e5-1e6) otherwise produce FPC scores ~1e11 whose variance-component fit goes
  singular, forcing genotype variance to 0 and a spurious H2~0; standardizing
  fixes this (`side_area` FPC1 H2 0.007 → 0.30) without changing well-conditioned
  traits. Diagnostic + eigenfunctions: `4_FPCA/script/4b_fpca_eigen_side_area.R`
  → `fpca_eigen_side_area.{png,csv}` (eigenfunctions on a 500-pt day grid).
  Remaining near-zero FPC H2 traits (ellipse_angle, center_of_mass_x, …) are
  genuinely non-heritable positional/orientation features, not artifacts.
- LDA LD1 leave-one-genotype-out CV heritability ≈ **0.666**.
