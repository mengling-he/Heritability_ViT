# =============================================================================
# 5b_lda_cv5_genotype.R
#   5-fold genotype-level CV LDA:
#   - split the 60 genotypes into 5 folds (12 genotypes each)
#   - per fold: scale with training stats, fit LDA on the 48 training genotypes,
#     project the 12 held-out genotypes -> out-of-fold LD1-5 scores
#   - H2 (date-as-covariate model) of the held-out LD scores, computed within
#     each fold on the 12 test genotypes; summarised over folds x repeats
#   - optional repeats with different random fold splits (n_repeats)
#
# Notes:
# - Folds split GENOTYPES (all replicates of a genotype stay in the same fold),
#   so the LDA never sees the genotypes it is scored on.
# - Per-fold H2 uses each fold's own LD axes (LDk = k-th discriminant of that
#   fold's model), so no cross-fold alignment is needed for H2.
# - Fold LDs ARE matched + sign-aligned to a full-data reference LDA, but only
#   to average loadings across folds (same procedure as 5a_lda_logo.R).
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(MASS)
  library(ggplot2)
  library(viridis)
  library(here)
})
select <- dplyr::select

source(here("renv", "activate.R"))
segmentation_method <- "vit"  # Change to "unet" when needed
base_dir <- switch(
  segmentation_method,
  "unet" = "Heritability_Unet",
  "vit"  = "Heritability_ViT",
  stop("Unknown segmentation method")
)
code_path <- function(...) here(base_dir, "Code", ...)
prep_path <- function(...) here(base_dir, "Analysis", "1_preprocessing", "result", ...)
h2_path <- function(...) here(base_dir, "Analysis", "2_heritability", "result", ...)
result_path <- function(...) here(base_dir, "Analysis", "5_LDA", "result", ...)
RESULT_DIR <- result_path()
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
source(code_path("heritability_fun.R"))

# ---- CV settings ----
n_folds   <- 5
# >1 = repeated 5-fold CV with different random genotype splits;
# override from the command line, e.g. Rscript 5b_lda_cv5_genotype.R 1
args      <- commandArgs(trailingOnly = TRUE)
n_repeats <- 10
#n_repeats <- if (length(args) >= 1) as.integer(args[1]) else 5
base_seed <- 123   # repeat r uses seed base_seed + r - 1
n_dims    <- 5
lda_dim_names <- paste0("LD", 1:n_dims)

# all CV outputs go to a subfolder named after the settings, e.g. result/cv5_rep5/
cv_path <- function(...) result_path(sprintf("cv%d_rep%d", n_folds, n_repeats), ...)
dir.create(cv_path(), showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Helper functions (same as 5a_lda_logo.R)
# =============================================================================
cosine_sim <- function(a, b) {
  denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (denom == 0) return(NA_real_)
  sum(a * b) / denom
}

# all permutations of integer vector; fine here because n_dims = 5 (120 perms)
permute_indices <- function(x) {
  if (length(x) == 1) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], permute_indices(x[-i]))
  }))
}

# best one-to-one matching from fold LDs to reference LDs (max total |cosine|)
best_ld_match <- function(ref_scaling, fold_scaling, perm_mat) {
  sim_mat <- abs(outer(seq_len(ncol(ref_scaling)), seq_len(ncol(fold_scaling)),
                       Vectorize(function(i, j) cosine_sim(ref_scaling[, i], fold_scaling[, j]))))
  dimnames(sim_mat) <- list(colnames(ref_scaling), colnames(fold_scaling))
  perm_scores <- apply(perm_mat, 1, function(p) sum(sim_mat[cbind(seq_len(nrow(sim_mat)), p)]))
  list(sim_mat = sim_mat, best_perm = perm_mat[which.max(perm_scores), ])
}

scale_with <- function(df, cols, means, sds) {
  df[cols] <- scale(df[cols], center = means, scale = sds)
  df
}

# =============================================================================
# Read and preprocess data (same as 5a_lda_logo.R)
# =============================================================================
features_df0_top_side <- read.csv(
  prep_path(sprintf("side_top_features_%s_linearimputed_day.csv", segmentation_method)),
  stringsAsFactors = FALSE
)
slow_growers <- strsplit(trimws(readLines(prep_path("slow_grower.txt"), warn = FALSE)), ",")[[1]]
df_cv <- features_df0_top_side %>% filter(!Plant.Info %in% slow_growers) %>%
  tidyr::extract(Plant.Info, into = c("Genotype", "Replicate"),
                 regex = "^\\d+_(.+)_(\\d+)$", remove = FALSE, convert = FALSE)
df_cv$date <- as.Date(df_cv$date)
df_cv$Day  <- as.numeric(df_cv$date - min(df_cv$date))
df_cv$Genotype <- factor(df_cv$Genotype)

traits <- grep("^(side|top)_", names(df_cv), value = TRUE)
nzv <- sapply(df_cv[traits], function(x) { v <- var(x, na.rm = TRUE); is.finite(v) && v > 1e-8 })
metric_cols_lda <- traits[nzv]

genotypes <- sort(unique(as.character(df_cv$Genotype)))
cat(sprintf("5-fold genotype CV LDA: %d traits, %d genotypes (~%d held out per fold), %d repeat(s)\n",
            length(metric_cols_lda), length(genotypes), round(length(genotypes) / n_folds), n_repeats))

# =============================================================================
# Reference LDA (full data) — used only to match / sign-align fold loadings
# =============================================================================
ref_means <- sapply(df_cv[metric_cols_lda], mean, na.rm = TRUE)
ref_sds   <- sapply(df_cv[metric_cols_lda], sd,   na.rm = TRUE)
ref_sds[ref_sds == 0] <- 1
lda_ref <- lda(Genotype ~ ., data = scale_with(df_cv, metric_cols_lda, ref_means, ref_sds)[, c("Genotype", metric_cols_lda)])
ref_scaling <- lda_ref$scaling[, lda_dim_names, drop = FALSE]
perm_mat    <- permute_indices(seq_along(lda_dim_names))

# =============================================================================
# Repeated 5-fold CV loop
# =============================================================================
fold_h2_rows   <- list()
loading_rows   <- list()
fold_map_rows  <- list()

for (r in seq_len(n_repeats)) {
  # balanced random split: shuffle genotypes, deal them into folds (12 / fold)
  set.seed(base_seed + r - 1)
  fold_map <- tibble(Genotype = sample(genotypes),
                     Fold     = rep_len(seq_len(n_folds), length(genotypes)))
  fold_map_rows[[r]] <- fold_map %>% mutate(Repeat = r)

  for (k in seq_len(n_folds)) {
    test_genos <- fold_map$Genotype[fold_map$Fold == k]
    message(sprintf("=== Repeat %d, fold %d/%d: %d test genotypes ===", r, k, n_folds, length(test_genos)))

    train_df <- df_cv %>% filter(!Genotype %in% test_genos) %>% droplevels()
    test_df  <- df_cv %>% filter(Genotype %in% test_genos)

    # scale using training statistics only
    train_means <- sapply(train_df[metric_cols_lda], mean, na.rm = TRUE)
    train_sds   <- sapply(train_df[metric_cols_lda], sd,   na.rm = TRUE)
    train_sds[train_sds == 0] <- 1
    train_scaled <- scale_with(train_df, metric_cols_lda, train_means, train_sds)
    test_scaled  <- scale_with(test_df,  metric_cols_lda, train_means, train_sds)

    # fit LDA on the 48 training genotypes, project the 12 held-out genotypes
    lda_model <- lda(Genotype ~ ., data = train_scaled[, c("Genotype", metric_cols_lda)])
    fold_scaling_raw <- lda_model$scaling[, lda_dim_names, drop = FALSE]
    lda_x_test_raw <- as.data.frame(predict(lda_model, newdata = test_scaled)$x)[, lda_dim_names, drop = FALSE]

    meta <- test_scaled[, c("Plant.Info", "date", "Day", "Genotype", "Replicate")]
    meta$Genotype <- droplevels(factor(meta$Genotype))

    # per-fold H2 on the fold's own LD axes (12 test genotypes)
    h2_k <- calc_H2_fixedtime(cbind(meta, lda_x_test_raw), lda_dim_names, day_col = "Day")
    fold_h2_rows[[length(fold_h2_rows) + 1]] <- h2_k %>% mutate(Repeat = r, Fold = k)

    # match fold LDs to reference LDs, then sign-align (for averaging loadings)
    match_obj    <- best_ld_match(ref_scaling, fold_scaling_raw, perm_mat)
    matched      <- colnames(fold_scaling_raw)[match_obj$best_perm]
    fold_scaling <- fold_scaling_raw[, matched, drop = FALSE]
    colnames(fold_scaling) <- lda_dim_names
    ld_signs <- sapply(lda_dim_names, function(d) {
      if (cosine_sim(ref_scaling[, d], fold_scaling[, d]) < 0) -1 else 1
    })
    fold_scaling <- sweep(fold_scaling, 2, ld_signs, `*`)

    loading_rows[[length(loading_rows) + 1]] <- as.data.frame(fold_scaling) %>%
      rownames_to_column("Variable") %>%
      mutate(Repeat = r, Fold = k,
             Matched_Fold_LD = paste(matched, collapse = ","),
             AbsCosine_LD1 = abs(cosine_sim(ref_scaling[, "LD1"], fold_scaling[, "LD1"])))
  }
}

# =============================================================================
# Collect and save results
# =============================================================================
fold_h2  <- bind_rows(fold_h2_rows) %>% mutate(Trait = factor(Trait, levels = lda_dim_names))
loadings <- bind_rows(loading_rows)

# summary over all folds x repeats, with a 95% CI for the mean per-fold H2 from the
# corrected repeated-CV t-interval (Nadeau & Bengio 2003; Bouckaert & Frank 2004):
#   SE = SD * sqrt(1/J + n_test/n_train),  CI = mean +/- t(0.975, J-1) * SE
# J = folds x repeats. The n_test/n_train term (12/48 genotypes) corrects for the
# overlapping training sets; a naive SD/sqrt(J) would be far too narrow.
test_train_ratio <- (length(genotypes) / n_folds) / (length(genotypes) - length(genotypes) / n_folds)
h2_summary <- fold_h2 %>%
  group_by(Trait) %>%
  summarise(H2_mean   = mean(H2,   na.rm = TRUE),
            H2_sd     = sd(H2,     na.rm = TRUE),
            H2_median = median(H2, na.rm = TRUE),
            H2_min    = min(H2,    na.rm = TRUE),
            H2_max    = max(H2,    na.rm = TRUE),
            n_folds   = sum(!is.na(H2)), .groups = "drop") %>%
  mutate(H2_se_corrected = H2_sd * sqrt(1 / n_folds + test_train_ratio),
         CI95_low = H2_mean - qt(0.975, df = n_folds - 1) * H2_se_corrected,
         CI95_up  = H2_mean + qt(0.975, df = n_folds - 1) * H2_se_corrected) %>%
  arrange(Trait)

# mean per-fold H2 within each repeat (how much the random split matters), with a
# 95% t-interval over that repeat's folds: mean +/- t(0.975, n-1) * SD / sqrt(n).
# Test genotypes are disjoint across folds, but training sets overlap, so the
# fold H2s are not fully independent and this CI is somewhat optimistic.
h2_by_repeat <- fold_h2 %>%
  group_by(Repeat, Trait) %>%
  summarise(n_folds = sum(!is.na(H2)),
            H2_mean = mean(H2, na.rm = TRUE),
            H2_sd   = sd(H2,   na.rm = TRUE), .groups = "drop") %>%
  mutate(H2_se    = H2_sd / sqrt(n_folds),
         CI95_low = H2_mean - qt(0.975, df = n_folds - 1) * H2_se,
         CI95_up  = H2_mean + qt(0.975, df = n_folds - 1) * H2_se) %>%
  arrange(Repeat, Trait)

write.csv(bind_rows(fold_map_rows), cv_path("cv5_LDA_fold_assignment.csv"), row.names = FALSE)
write.csv(fold_h2,    cv_path(sprintf("cv5_LDA_h2_per_fold_%s.csv", segmentation_method)), row.names = FALSE)
write.csv(h2_summary, cv_path(sprintf("cv5_LDA_h2_summary_%s.csv",  segmentation_method)), row.names = FALSE)
write.csv(h2_by_repeat, cv_path(sprintf("cv5_LDA_h2_by_repeat_%s.csv", segmentation_method)), row.names = FALSE)

cat("\n5-fold genotype CV LDA heritability:\n")
print(as.data.frame(h2_summary), row.names = FALSE, digits = 3)

# =============================================================================
# Plots
# =============================================================================
# H2 per LD: mean per-fold H2 with corrected 95% CI (see h2_summary)
h2_summary <- read.csv(cv_path(sprintf("cv5_LDA_h2_summary_%s.csv", segmentation_method)))
p_h2 <- ggplot(h2_summary, aes(Trait, H2_mean)) +
  geom_col(fill = "salmon", width = 0.65) +
  geom_errorbar(aes(ymin = CI95_low, ymax = CI95_up), width = 0.25) +
  geom_text(aes(y = 0, label = sprintf("%.3f", H2_mean)), vjust = -0.6, size = 3.5) +
  coord_cartesian(ylim = c(0, 0.8)) +   # clip (not drop) CI bars that cross 0 or 1
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = expression(H^2))
ggsave(cv_path(sprintf("cv5_LDA_h2_%s2.png", segmentation_method)), p_h2, width = 4, height = 6, dpi = 300)

# same bar plot without error bars
p_h2_noci <- ggplot(h2_summary, aes(Trait, H2_mean)) +
  geom_col(fill = "salmon", width = 0.65) +
  geom_text(aes(y = 0, label = sprintf("%.3f", H2_mean)), vjust = -0.6, size = 3.5) +
  coord_cartesian(ylim = c(0, 0.7)) +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = expression(H^2))
ggsave(cv_path(sprintf("cv5_LDA_h2_noCI_%s_@.png", segmentation_method)), p_h2_noci, width = 7, height = 5, dpi = 300)

# distribution of per-fold H2: boxplot over all folds, dots colored by repeat
p_folds <- ggplot(fold_h2, aes(Trait, H2)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey95") +
  geom_point(aes(color = factor(Repeat)), size = 1.8, alpha = 0.7,
             position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = expression(H^2), color = "Repeat")
ggsave(cv_path(sprintf("cv5_LDA_h2_folds_%s.png", segmentation_method)), p_folds, width = 7, height = 4.5, dpi = 300)

# top-20 mean aligned loadings, bar fill = raw trait H2 (same style as 5a)
raw_h2 <- read.csv(h2_path(sprintf("h2_timefixed_day_%s.csv", segmentation_method)),
                   stringsAsFactors = FALSE) %>% select(Trait, H2)
mean_loadings <- loadings %>%
  group_by(Variable) %>%
  summarise(across(all_of(lda_dim_names), \(x) mean(x, na.rm = TRUE), .names = "mean_{.col}"),
            across(all_of(lda_dim_names), \(x) sd(x,   na.rm = TRUE), .names = "sd_{.col}"),
            .groups = "drop") %>%
  left_join(raw_h2, by = c("Variable" = "Trait"))
write.csv(mean_loadings, cv_path("cv5_LDA_mean_loadings.csv"), row.names = FALSE)

make_loadplot <- function(ld, top_n = 20) {
  col_mean <- paste0("mean_", ld)
  top <- mean_loadings %>% arrange(desc(abs(.data[[col_mean]]))) %>% slice_head(n = top_n)
  ggplot(top, aes(reorder(Variable, .data[[col_mean]]), .data[[col_mean]], fill = H2)) +
    geom_col() + coord_flip() +
    scale_fill_viridis(option = "plasma", name = expression(H^2)) +
    theme_minimal(base_size = 16) +
    theme(axis.text.y = element_text(size = 14)) +
    labs(x = NULL, y = sprintf("Mean loading across CV folds (%s)", ld))
}
for (ld in c("LD1", "LD2")) {
  ggsave(cv_path(sprintf("cv5_LDA_loadings_%s_h2color_%s.png", ld, segmentation_method)),
         make_loadplot(ld), width = 7, height = 5.5, dpi = 300)
}

cat("\n5-fold genotype CV LDA finished. Outputs written to:\n", cv_path(), "\n")
