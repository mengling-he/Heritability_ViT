suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(MASS); library(lme4); library(ggplot2); library(viridis)
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

# =============================================================================
# Data: daily-aligned, slow growers removed, complete cases, scaled traits
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
# trait set: drop near-constant
# metric_cols <- setdiff(colnames(features_df_top_side),
#                        c("Plant.Info", "date", "Genotype", "Replicate", "Day"))
# metric_cols_lda <- setdiff(metric_cols, c("top_blue_min", "top_saturation_max"))
# trait set: drop near-constant
traits <- grep("^(side|top)_", names(features_df_top_side), value = TRUE)
nzv <- sapply(features_df_top_side[traits], function(x) { v <- var(x, na.rm = TRUE); is.finite(v) && v > 1e-8 })
traits <- traits[nzv]
metric_cols_lda <- traits
# exclude two genotypes (same as original script)
# df_logo <- features_df_top_side %>%
#   filter(!Genotype %in% c("BESC-1018", "KLNA-20-3"))

# origin_date <- min(features_df_top_side$date)
n_dims      <- 5   # first 5 LDA directions
lda_dim_names <- paste0("LD", 1:n_dims)

df_logo <- features_df_top_side
# =========================================================
# Reference LDA (full data) — used only for sign alignment
# Sign flips across LOGO folds would destroy between-genotype
# variance when pooling projections, so we anchor all folds
# to a common reference direction.
# =========================================================
ref_means <- sapply(df_logo[metric_cols_lda], mean, na.rm = TRUE)
ref_sds   <- sapply(df_logo[metric_cols_lda], sd,   na.rm = TRUE)
ref_sds[ref_sds == 0] <- 1

df_logo_ref <- df_logo
df_logo_ref[metric_cols_lda] <- scale(df_logo[metric_cols_lda], center = ref_means, scale = ref_sds)

# keep <- sapply(df_logo_ref[, metric_cols_lda], function(x) var(x, na.rm = TRUE) > 0)
# metric_cols_lda_clean <- metric_cols_lda[keep]

lda_ref <- lda(Genotype ~ ., data = df_logo_ref[, c("Genotype", metric_cols_lda)])
ref_scaling    <- lda_ref$scaling[, lda_dim_names, drop = FALSE]  # reference loading vectors

# =========================================================
# LOGO loop
# =========================================================
genotypes <- sort(unique(df_logo$Genotype))
n_geno    <- length(genotypes)

pooled_projections <- list()
all_loadings       <- list()
all_sc             <- list()

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

  # sign-align each LD to the reference via dot product of loading vectors
  fold_scaling <- lda_model$scaling[, lda_dim_names, drop = FALSE]
  ld_signs <- sapply(lda_dim_names, function(d) {
    if (sum(ref_scaling[, d] * fold_scaling[, d]) < 0) -1L else 1L
  })

  # project left-out genotype onto fold LDA space; apply sign alignment
  lda_x_test <- data.frame(predict(lda_model, newdata = test_scaled)$x)[, lda_dim_names, drop = FALSE]
  for (d in lda_dim_names) {
    lda_x_test[[d]] <- lda_x_test[[d]] * ld_signs[[d]]
  }

  proj_df <- cbind(
    test_scaled[, c("date", "Day","Genotype", "Replicate", "Plant.Info")],
    lda_x_test
  )
  proj_df$date <- as.Date(proj_df$date)
  proj_df$Day  <- as.numeric(proj_df$Day)
  pooled_projections[[i]] <- proj_df

  # loadings (sign-aligned)
  ld_df <- data.frame(sweep(fold_scaling, 2, ld_signs, `*`))
  ld_df$Variable <- rownames(ld_df)
  ld_df$Left_out <- g_out
  all_loadings[[i]] <- ld_df

  # structure coefficients on training set (sign-aligned)
  lda_x_train <- data.frame(predict(lda_model, newdata = train_scaled)$x)[, lda_dim_names, drop = FALSE]
  for (d in lda_dim_names) {
    lda_x_train[[d]] <- lda_x_train[[d]] * ld_signs[[d]]
  }
  sc_mat <- cor(train_scaled[, metric_cols_lda], lda_x_train)
  sc_df  <- data.frame(sc_mat)
  sc_df$Variable <- rownames(sc_df)
  sc_df$Left_out <- g_out
  all_sc[[i]] <- sc_df
}

# =========================================================
# Pool all left-out projections and compute H²
# =========================================================
pooled_df <- bind_rows(pooled_projections)

h2_logo <- calc_H2_fixedtime(pooled_df, lda_dim_names, day_col = "Day")

# =========================================================
# Aggregate loadings and structure coefficients across folds
# =========================================================
loadings_all <- bind_rows(all_loadings)

loadings_summary_LD1 <- loadings_all %>%
  group_by(Variable) %>%
  summarise(mean_LD1     = mean(LD1, na.rm = TRUE),
            mean_abs_LD1 = mean(abs(LD1), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_LD1)) %>%
  slice(1:20)

loadings_summary_LD2 <- loadings_all %>%
  group_by(Variable) %>%
  summarise(mean_LD2     = mean(LD2, na.rm = TRUE),
            mean_abs_LD2 = mean(abs(LD2), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_LD2)) %>%
  slice(1:20)

sc_all <- bind_rows(all_sc)

sc_summary_LD1 <- sc_all %>%
  group_by(Variable) %>%
  summarise(mean_SC_LD1 = mean(LD1, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(abs(mean_SC_LD1))) %>%
  slice(1:20)

# =========================================================
# Plots
# =========================================================

p_h2_logo <- ggplot(h2_logo, aes(x = reorder(Trait, H2), y = H2)) +
  geom_col(fill = "salmon") +
  geom_text(aes(label = round(H2, 3)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  theme_minimal() +
  labs(title = "LOGO CV LDA: H² on pooled left-out projections (first 5 LDs)",
       x = "", y = "H²")

p_ld1_logo <- ggplot(loadings_summary_LD1, aes(x = reorder(Variable, mean_LD1), y = mean_LD1)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "LOGO CV LDA: Mean Loadings for LD1 (top 20)",
       x = "", y = "Mean loading across LOGO folds")

p_ld2_logo <- ggplot(loadings_summary_LD2, aes(x = reorder(Variable, mean_LD2), y = mean_LD2)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "LOGO CV LDA: Mean Loadings for LD2 (top 20)",
       x = "", y = "Mean loading across LOGO folds")

h2_raw <- read.csv(h2_path(sprintf("h2_timefixed_day_%s.csv", segmentation_method)))

sc_ld1_h2 <- sc_summary_LD1 %>%
  left_join(h2_raw %>% dplyr::select(Trait, H2), by = c("Variable" = "Trait"))

p_sc_ld1_logo <- ggplot(sc_ld1_h2, aes(x = reorder(Variable, abs(mean_SC_LD1)),
                                         y = mean_SC_LD1, fill = H2)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "LOGO CV LDA: Mean Structure Coefficients (LD1, top 20)",
       x = "", y = "Mean SC with LD1", fill = "H²")

# =========================================================
# Save outputs
# =========================================================
write.csv(h2_logo,   result_path("logo_LDA_h2.csv"),                 row.names = FALSE)
write.csv(pooled_df, result_path("logo_LDA_pooled_projections.csv"),  row.names = FALSE)

ggsave(result_path("logo_LDA_h2.png"),      p_h2_logo,    width = 7, height = 5, dpi = 300)
ggsave(result_path("logo_LDA_LD1.png"),     p_ld1_logo,   width = 7, height = 5, dpi = 300)
ggsave(result_path("logo_LDA_LD2.png"),     p_ld2_logo,   width = 7, height = 5, dpi = 300)
ggsave(result_path("logo_LDA_sc_ld1.png"),  p_sc_ld1_logo, width = 7, height = 5, dpi = 300)

message("LOGO CV complete. Files saved to ", RESULT_DIR)
