# Heritability_ViT

Broad-sense heritability (H²) of longitudinal poplar image phenotypes, using
side- and top-view features extracted from **ViT segmentation** and compared
against a **thresholding** baseline.

## Folder structure

```
raw_data/       ViT and thresholding side/top feature CSVs
Code/           shared functions (heritability_fun.R, plot_fun.R, FPCA_code/)
Analysis/
  1_preprocessing/         slow growers, side+top merge, daily imputation
  2_heritability/          per-date and date-as-covariate H²
  3_threshold_comparison/  ViT vs thresholding H²
  4_FPCA/                  functional PCA scores and their H²
  5_LDA/                   LDA scores and their H² (LOGO and 5-fold CV)
```

Each step has a `script/` folder; outputs are written to `result/`.

## Workflow

| Step | Script | Purpose |
|---|---|---|
| 1 | `1_slow_grower.R` | Shared slow-grower list (side-area growth < 5%) |
| 1 | `2_preprocess.R` / `2a_preprocess_thresholding.R` | Merge views, impute onto a daily grid (ViT / thresholding) |
| 2 | `2_heritability.R` / `2a_heritabiliyy_threshold.R` | H² of every trait (ViT / thresholding) |
| 2 | `2b_view_figures.R` | Side vs top H² figures |
| 3 | `3_compare_threshold.R` | ViT vs thresholding H² comparison |
| 4 | `4_fpca.R` | FPC scores and their H² |
| 5 | `5_lda.R` | LDA on all genotypes |
| 5 | `5a_lda_logo.R` | Leave-one-genotype-out LDA, H² with bootstrap CI |
| 5 | `5b_lda_cv5_genotype.R` | Repeated 5-fold genotype CV LDA, per-fold H² with CI |

Run the steps in order; later steps read earlier results.

## Usage

Scripts use `here()` and the renv library of the parent `Heritability/`
project, so run them from that project root:

```bash
Rscript Heritability_ViT/Analysis/1_preprocessing/script/1_slow_grower.R
Rscript Heritability_ViT/Analysis/1_preprocessing/script/2_preprocess.R
# ... then steps 2 to 5
Rscript Heritability_ViT/Analysis/5_LDA/script/5b_lda_cv5_genotype.R 10   # 10 repeats
```

## Models

| Analysis | Model |
|---|---|
| Per-date H² | `trait ~ (1 \| Genotype)` |
| Date-as-covariate H² | `trait ~ Day + (1 \| Genotype) + (1 \| Genotype:Replicate)` |
| FPC score H² | `score ~ (1 \| Genotype)` |

H² = σ²_G / (σ²_G + σ²_E). Genotype and replicate are parsed from `Plant.Info`.
