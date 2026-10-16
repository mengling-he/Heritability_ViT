# =============================================================================
#   1) LOGO-CV LDA with fold-to-reference dimension matching
#   2) sign alignment and axis-switching diagnostics
#   3) loading stability summaries / plots
#   4) H2 of pooled LOGO LD1-5 scores
#   6) H2-colored mean loading plots
#
# Notes:
# - This script is intended to replace the current LDA_LOGO_loadings_h2color.R
# - The bootstrap CI implemented here is a fast/practical CI that conditions on
#   the already learned pooled LOGO projections. It does NOT rerun the entire
#   LDA + LOGO pipeline inside each bootstrap replicate.
# - Select the traits that are not zero variance, include all genotypes
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(MASS)
  library(ggplot2)
  library(viridis)
  library(here)
  library(purrr)
  library(stringr)
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

# =============================================================================
# Helper functions
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
    rest <- x[-i]
    p <- permute_indices(rest)
    cbind(x[i], p)
  }))
}

loading_topk_overlap <- function(v1, v2, k = 20) {
  k <- min(k, length(v1), length(v2))
  s1 <- names(sort(abs(v1), decreasing = TRUE))[seq_len(k)]
  s2 <- names(sort(abs(v2), decreasing = TRUE))[seq_len(k)]
  length(intersect(s1, s2)) / k
}

# find the best one-to-one matching from fold LDs to reference LDs
best_ld_match <- function(ref_scaling, fold_scaling, perm_mat) {
  ref_names  <- colnames(ref_scaling)
  fold_names <- colnames(fold_scaling)

  sim_mat <- matrix(NA_real_, nrow = length(ref_names), ncol = length(fold_names),
                    dimnames = list(ref_names, fold_names))

  for (i in seq_along(ref_names)) {
    for (j in seq_along(fold_names)) {
      sim_mat[i, j] <- abs(cosine_sim(ref_scaling[, ref_names[i]],
                                      fold_scaling[, fold_names[j]]))
    }
  }

  perm_scores <- apply(perm_mat, 1, function(p) {
    sum(sim_mat[cbind(seq_len(nrow(sim_mat)), p)])
  })

  best_idx  <- which.max(perm_scores)
  best_perm <- perm_mat[best_idx, ]

  list(
    sim_mat    = sim_mat,
    best_perm  = best_perm,
    best_score = perm_scores[best_idx]
  )
}

# =============================================================================
# Read and preprocess data
# =============================================================================
features_df0_top_side <- read.csv(
  prep_path(sprintf("side_top_features_%s_linearimputed_day.csv",segmentation_method)),
  stringsAsFactors = FALSE
)
slow_growers <- strsplit(trimws(readLines(
  prep_path("slow_grower.txt"),
  warn = FALSE
)),
","
)[[1]]
features_df_top_side <- features_df0_top_side %>% filter(!Plant.Info %in% slow_growers) %>%
  tidyr::extract(Plant.Info, into = c("Genotype", "Replicate"),
                 regex = "^\\d+_(.+)_(\\d+)$", remove = FALSE, convert = FALSE)
features_df_top_side$date <- as.Date(features_df_top_side$date)
features_df_top_side$Day  <- as.numeric(features_df_top_side$date - min(features_df_top_side$date))
features_df_top_side$Genotype <- factor(features_df_top_side$Genotype)

traits <- grep("^(side|top)_", names(features_df_top_side), value = TRUE)
nzv <- sapply(features_df_top_side[traits], function(x) { v <- var(x, na.rm = TRUE); is.finite(v) && v > 1e-8 })
traits <- traits[nzv]
metric_cols_lda <- traits

df_logo <- features_df_top_side
origin_date   <- min(df_logo$date)
n_dims        <- 5
lda_dim_names <- paste0("LD", 1:n_dims)
perm_mat      <- permute_indices(seq_along(lda_dim_names))

# =============================================================================
# Reference LDA (full data) — used for fold matching and sign alignment
# =============================================================================
ref_means <- sapply(df_logo[metric_cols_lda], mean, na.rm = TRUE)
ref_sds   <- sapply(df_logo[metric_cols_lda], sd,   na.rm = TRUE)
ref_sds[ref_sds == 0] <- 1

df_logo_ref <- df_logo
df_logo_ref[metric_cols_lda] <- scale(df_logo[metric_cols_lda], center = ref_means, scale = ref_sds)

lda_ref <- lda(Genotype ~ ., data = df_logo_ref[, c("Genotype", metric_cols_lda)])
ref_scaling <- lda_ref$scaling[, lda_dim_names, drop = FALSE]

write.csv(
  cbind(Variable = rownames(ref_scaling), as.data.frame(ref_scaling)),
  result_path("logo_LDA_reference_loadings.csv"),
  row.names = FALSE
)

# =============================================================================
# LOGO loop — fit LDA leaving out each genotype in turn
# =============================================================================

genotypes <- sort(unique(df_logo$Genotype))
n_geno    <- length(genotypes)

all_loadings        <- list()
pooled_projections  <- list()
stability_rows      <- list()
sim_mats_long       <- list()
matching_rows       <- list()

for (i in seq_along(genotypes)) {
  g_out <- genotypes[i]
  message("=== LOGO fold ", i, "/", n_geno, ": leaving out ", g_out, " ===")

  train_df <- df_logo %>% filter(Genotype != g_out)
  test_df  <- df_logo %>% filter(Genotype == g_out)

  # scale using training statistics only
  train_means <- sapply(train_df[metric_cols_lda], mean, na.rm = TRUE)
  train_sds   <- sapply(train_df[metric_cols_lda], sd,   na.rm = TRUE)
  train_sds[train_sds == 0] <- 1

  train_scaled <- train_df
  test_scaled  <- test_df
  train_scaled[metric_cols_lda] <- scale(train_df[metric_cols_lda], center = train_means, scale = train_sds)
  test_scaled[metric_cols_lda]  <- scale(test_df[metric_cols_lda],  center = train_means, scale = train_sds)

  # fit LDA without g_out
  lda_model <- lda(Genotype ~ ., data = train_scaled[, c("Genotype", metric_cols_lda)])
  fold_scaling_raw <- lda_model$scaling[, lda_dim_names, drop = FALSE]

  # best one-to-one matching of fold LDs to reference LDs
  match_obj <- best_ld_match(ref_scaling, fold_scaling_raw, perm_mat)
  matched_fold_names <- colnames(fold_scaling_raw)[match_obj$best_perm]

  matching_rows[[i]] <- tibble(
    Left_out        = g_out,
    Ref_LD          = lda_dim_names,
    Matched_Fold_LD = matched_fold_names,
    Switched        = matched_fold_names != lda_dim_names
  )

  # reorder fold LDs to align with reference LD1...LD5
  fold_scaling <- fold_scaling_raw[, matched_fold_names, drop = FALSE]
  colnames(fold_scaling) <- lda_dim_names

  # sign-align AFTER matching
  ld_signs <- sapply(lda_dim_names, function(d) {
    if (cosine_sim(ref_scaling[, d], fold_scaling[, d]) < 0) -1L else 1L
  })
  fold_scaling_aligned <- sweep(fold_scaling, 2, ld_signs, `*`)

  # save aligned loadings
  ld_df <- as.data.frame(fold_scaling_aligned)
  ld_df$Variable <- rownames(fold_scaling_aligned)
  ld_df$Left_out <- g_out
  all_loadings[[i]] <- ld_df

  # store similarity matrix in long form for diagnostics
  sim_df <- as.data.frame(as.table(match_obj$sim_mat), stringsAsFactors = FALSE)
  colnames(sim_df) <- c("Ref_LD", "Fold_LD", "AbsCosine")
  sim_df$Left_out <- g_out
  sim_mats_long[[i]] <- sim_df

  # per-fold stability summary after matching + sign alignment
  stab_i <- bind_rows(lapply(seq_along(lda_dim_names), function(k) {
    ref_d <- lda_dim_names[k]
    matched_orig <- matched_fold_names[k]
    tibble(
      Left_out        = g_out,
      Ref_LD          = ref_d,
      Matched_Fold_LD = matched_orig,
      SignedCosine    = cosine_sim(ref_scaling[, ref_d], fold_scaling_aligned[, ref_d]),
      AbsCosine       = abs(cosine_sim(ref_scaling[, ref_d], fold_scaling_aligned[, ref_d])),
      Top20Overlap    = loading_topk_overlap(ref_scaling[, ref_d], fold_scaling_aligned[, ref_d], k = 20),
      Switched        = matched_orig != ref_d
    )
  }))
  stability_rows[[i]] <- stab_i

  # project left-out genotype onto fold LDA space and apply same matching/sign alignment
  lda_x_test_raw <- data.frame(predict(lda_model, newdata = test_scaled)$x)[, colnames(fold_scaling_raw), drop = FALSE]
  lda_x_test <- lda_x_test_raw[, matched_fold_names, drop = FALSE]
  colnames(lda_x_test) <- lda_dim_names
  for (d in lda_dim_names) {
    lda_x_test[[d]] <- lda_x_test[[d]] * ld_signs[[d]]
  }

  proj_df <- cbind(
    test_scaled[, c("date", "Genotype", "Replicate", "Plant.Info")],
    lda_x_test
  )
  proj_df$date <- as.Date(proj_df$date)
  proj_df$Day  <- as.numeric(proj_df$date - origin_date)
  pooled_projections[[i]] <- proj_df
}

# =============================================================================
# Stability diagnostics
# =============================================================================

loadings_all   <- bind_rows(all_loadings)
stability_all  <- bind_rows(stability_rows)
sim_all        <- bind_rows(sim_mats_long)
matching_all   <- bind_rows(matching_rows)

write.csv(stability_all, result_path("logo_LDA_loading_stability.csv"), row.names = FALSE)
write.csv(sim_all,       result_path("logo_LDA_similarity_matrices_long.csv"), row.names = FALSE)
write.csv(matching_all,  result_path("logo_LDA_matching_by_fold.csv"), row.names = FALSE)

stability_summary <- stability_all %>%
  group_by(Ref_LD) %>%
  summarise(
    MeanSignedCosine = mean(SignedCosine, na.rm = TRUE),
    MeanAbsCosine    = mean(AbsCosine, na.rm = TRUE),
    MinAbsCosine     = min(AbsCosine, na.rm = TRUE),
    MedianAbsCosine  = median(AbsCosine, na.rm = TRUE),
    MeanTop20Overlap = mean(Top20Overlap, na.rm = TRUE),
    SwitchRate       = mean(Switched, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(stability_summary, result_path("logo_LDA_loading_stability_summary.csv"), row.names = FALSE)

# plot: abs cosine stability by LD
p_stability_cos <- ggplot(stability_all, aes(x = Ref_LD, y = AbsCosine, fill = Ref_LD)) +
  geom_boxplot(alpha = 0.8, width = 0.65, outlier.size = 1) +
  geom_hline(yintercept = 0.9, linetype = 2, color = "gray40") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "|Cosine similarity to reference loading|")
ggsave(result_path("logo_LDA_stability_abs_cosine_boxplot.png"), p_stability_cos,
       width = 6.2, height = 4.6, dpi = 300)

# plot: top-20 overlap by LD
p_stability_top20 <- ggplot(stability_all, aes(x = Ref_LD, y = Top20Overlap, fill = Ref_LD)) +
  geom_boxplot(alpha = 0.8, width = 0.65, outlier.size = 1) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Top-20 absolute loading overlap")
ggsave(result_path("logo_LDA_stability_top20_overlap_boxplot.png"), p_stability_top20,
       width = 6.2, height = 4.6, dpi = 300)

# plot: switch rates
p_switch <- ggplot(stability_summary, aes(x = Ref_LD, y = SwitchRate, fill = Ref_LD)) +
  geom_col(width = 0.7) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Axis switching rate")
ggsave(result_path("logo_LDA_switch_rate.png"), p_switch,
       width = 6.0, height = 4.2, dpi = 300)

# average similarity matrix across folds
sim_mean <- sim_all %>%
  group_by(Ref_LD, Fold_LD) %>%
  summarise(MeanAbsCosine = mean(AbsCosine, na.rm = TRUE), .groups = "drop")

write.csv(sim_mean, result_path("logo_LDA_similarity_matrix_mean.csv"), row.names = FALSE)

p_sim_heat <- ggplot(sim_mean, aes(x = Fold_LD, y = Ref_LD, fill = MeanAbsCosine)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", MeanAbsCosine)), size = 3.5) +
  scale_fill_viridis(option = "magma", limits = c(0, 1), name = "Mean\n|cosine|") +
  theme_minimal(base_size = 11) +
  labs(x = "Fold LD", y = "Reference LD")
ggsave(result_path("logo_LDA_similarity_matrix_heatmap.png"), p_sim_heat,
       width = 6.0, height = 4.8, dpi = 300)

# =============================================================================
# Mean loadings across matched + sign-aligned LOGO folds
# =============================================================================

mean_loadings <- loadings_all %>%
  group_by(Variable) %>%
  summarise(
    across(all_of(lda_dim_names), \(x) mean(x, na.rm = TRUE), .names = "mean_{.col}"),
    across(all_of(lda_dim_names), \(x) sd(x, na.rm = TRUE),   .names = "sd_{.col}"),
    .groups = "drop"
  )

write.csv(mean_loadings, result_path("logo_LDA_mean_loadings.csv"), row.names = FALSE)

# =============================================================================
# H2 of LD1-5 scores (pooled left-out projections)
# =============================================================================

pooled_df <- bind_rows(pooled_projections)
write.csv(pooled_df, result_path("logo_LDA_pooled_projections.csv"), row.names = FALSE)

#pooled_df <- read.csv(result_path("logo_LDA_pooled_projections.csv"))
h2_logo <- calc_H2_fixedtime(pooled_df, lda_dim_names, day_col = "Day")
write.csv(h2_logo, result_path("logo_LDA_h2.csv"), row.names = FALSE)

# =============================================================================
# 95% CI for LOGO H2: genotype-level (cluster) bootstrap
#   resample the genotypes with replacement (all plants / days of a genotype
#   kept together), relabel duplicates as distinct genotypes, refit the H2 model.
#   Conditions on the LOGO projections: the 60 LDA models are NOT refit per
#   replicate, so LDA estimation uncertainty is not included.
# =============================================================================
B_boot    <- 200
boot_seed <- 2024
set.seed(boot_seed)
geno_ids  <- unique(as.character(pooled_df$Genotype))
by_geno   <- split(pooled_df, as.character(pooled_df$Genotype))

boot_h2 <- bind_rows(lapply(seq_len(B_boot), function(b) {
  if (b %% 20 == 0) message("LOGO bootstrap ", b, "/", B_boot)
  pick <- sample(geno_ids, length(geno_ids), replace = TRUE)
  boot_df <- bind_rows(lapply(seq_along(pick), function(i) {
    by_geno[[pick[i]]] %>% mutate(Genotype = paste0(pick[i], "_b", i))
  }))
  suppressMessages(suppressWarnings(
    calc_H2_fixedtime(boot_df, lda_dim_names, day_col = "Day")
  )) %>% select(Trait, H2) %>% mutate(Boot = b)
}))
write.csv(boot_h2, result_path("logo_LDA_h2_bootstrap_draws.csv"), row.names = FALSE)

h2_logo_ci <- h2_logo %>%
  select(Trait, H2) %>%
  left_join(boot_h2 %>% group_by(Trait) %>%
              summarise(H2_boot_mean = mean(H2, na.rm = TRUE),
                        H2_boot_se   = sd(H2,   na.rm = TRUE),
                        CI95_low     = unname(quantile(H2, 0.025, na.rm = TRUE)),
                        CI95_up      = unname(quantile(H2, 0.975, na.rm = TRUE)),
                        n_boot       = sum(!is.na(H2)), .groups = "drop"),
            by = "Trait")
write.csv(h2_logo_ci, result_path("logo_LDA_h2_CI.csv"), row.names = FALSE)
cat("\nLOGO LDA H2 with 95% genotype-bootstrap CI:\n")
print(as.data.frame(h2_logo_ci), row.names = FALSE, digits = 3)

p_h2_logo <- ggplot(h2_logo_ci, aes(x = reorder(Trait, -H2), y = H2)) +
  geom_col(fill = "salmon") +
  geom_errorbar(aes(ymin = CI95_low, ymax = CI95_up), width = 0.25) +
  geom_text(aes(y = 0, label = sprintf("%.3f", H2)), vjust = -0.6, size = 3.5) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  theme(plot.caption = element_text(size = 8, hjust = 0)) +
  labs(x = NULL, y = expression(H^2),
       caption = sprintf("Error bars: 95%% CI, genotype bootstrap (B = %d)", B_boot))
ggsave(result_path("logo_LDA_h2.png"), p_h2_logo, width = 4, height = 5, dpi = 300)


# =============================================================================
# Plots: top-20 mean LD1-LD4 loadings, bar fill = raw trait H2
# =============================================================================
raw_h2 <- read.csv(
  h2_path(sprintf("h2_timefixed_day_%s.csv",segmentation_method)),
  stringsAsFactors = FALSE
)%>%select(Trait, H2)

load_df <- mean_loadings %>%
  left_join(raw_h2, by = c("Variable" = "Trait"))

make_loadplot <- function(ld, top_n = 20) {
  col_mean <- paste0("mean_", ld)
  top <- load_df %>%
    arrange(desc(abs(.data[[col_mean]]))) %>%
    slice_head(n = top_n)

  ggplot(top, aes(reorder(Variable, .data[[col_mean]]), .data[[col_mean]], fill = H2)) +
    geom_col() +
    coord_flip() +
    scale_fill_viridis(option = "plasma", name = expression(H^2)) +
    theme_minimal(base_size = 16) +
    theme(axis.text.y = element_text(size = 14)) +
    labs(x = NULL, y = sprintf("Mean loading across LOGO folds (%s)", ld))
}

ggsave(result_path("logo_LDA_loadings_LD1_h2color.png"), make_loadplot("LD1"), width = 7, height = 5.5, dpi = 300)
ggsave(result_path("logo_LDA_loadings_LD2_h2color.png"), make_loadplot("LD2"), width = 7, height = 5.5, dpi = 300)
ggsave(result_path("logo_LDA_loadings_LD3_h2color.png"), make_loadplot("LD3"), width = 7, height = 5.5, dpi = 300)
ggsave(result_path("logo_LDA_loadings_LD4_h2color.png"), make_loadplot("LD4"), width = 7, height = 5.5, dpi = 300)
ggsave(result_path("logo_LDA_loadings_LD5_h2color.png"), make_loadplot("LD5"), width = 7, height = 5.5, dpi = 300)

# =============================================================================
# Optional: top loading tables
# =============================================================================

save_top_loadings <- function(ld, top_n = 20) {
  col_mean <- paste0("mean_", ld)
  out <- load_df %>%
    arrange(desc(abs(.data[[col_mean]]))) %>%
    slice_head(n = top_n) %>%
    select(Variable, all_of(col_mean), H2)
  write.csv(out, result_path(paste0("logo_LDA_top", top_n, "_", ld, "_loadings.csv")), row.names = FALSE)
}

invisible(lapply(lda_dim_names, save_top_loadings, top_n = 20))

cat("Complete LOGO-CV LDA analysis finished. Outputs written to:\n", RESULT_DIR, "\n")
