# =============================================================================
# LDA_LOGO_loadings_h2color.R
# Mean LOGO-CV LDA loadings (top 20 for LD1 / LD2), bar fill = raw trait H2,
# plus H2 of the LD1-5 scores themselves (pooled left-out projections).
# Combines the LOGO CV loop / sign-alignment logic from
# LDA_time_fixed_more_LOGO.R with the H2-colored loading plot style
# (make_loadplot) from 5_lda.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(MASS); library(ggplot2); library(viridis)
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

code_path   <- function(...) here(base_dir, "Code", ...)
prep_path   <- function(...) here(base_dir, "Analysis", "1_preprocessing", "result", ...)
h2_path     <- function(...) here(base_dir, "Analysis", "2_heritability", "result", ...)
result_path <- function(...) here(base_dir, "Analysis", "5_LDA", "result", ...)

RESULT_DIR <- result_path()
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
source(code_path("heritability_fun.R"))

# ---- read data ----
features_df0_top_side <- read.csv(prep_path(sprintf("side_top_features_%s_linearimputed_day.csv", segmentation_method)),
                                  stringsAsFactors = FALSE)

slow_growers <- strsplit(trimws(readLines(prep_path("slow_grower.txt"), warn = FALSE)), ",")[[1]]

features_df_top_side <- features_df0_top_side %>% filter(!Plant.Info %in% slow_growers) %>%
  tidyr::extract(Plant.Info, into = c("Genotype", "Replicate"),
                 regex = "^\\d+_(.+)_(\\d+)$", remove = FALSE, convert = FALSE)
features_df_top_side$date <- as.Date(features_df_top_side$date)
features_df_top_side$Day  <- as.numeric(features_df_top_side$date - min(features_df_top_side$date))

metric_cols <- setdiff(colnames(features_df_top_side),
                       c("Plant.Info", "date", "Genotype", "Replicate", "Day"))
metric_cols_lda <- setdiff(metric_cols, c("top_blue_min", "top_saturation_max"))

# exclude two genotypes (same as original LOGO script)
df_logo <- features_df_top_side %>%
  filter(!Genotype %in% c("BESC-1018", "KLNA-20-3"))

origin_date   <- min(df_logo$date)
n_dims        <- 5   # first 5 LDA directions
lda_dim_names <- paste0("LD", 1:n_dims)

# =========================================================
# Reference LDA (full data) — used only for sign alignment
# Sign flips across LOGO folds would cancel out when averaging
# loadings across folds, so anchor all folds to a common
# reference direction before aggregating.
# =========================================================
ref_means <- sapply(df_logo[metric_cols], mean, na.rm = TRUE)
ref_sds   <- sapply(df_logo[metric_cols], sd,   na.rm = TRUE)
ref_sds[ref_sds == 0] <- 1

df_logo_ref <- df_logo
df_logo_ref[metric_cols] <- scale(df_logo[metric_cols], center = ref_means, scale = ref_sds)

keep <- sapply(df_logo_ref[, metric_cols_lda], function(x) var(x, na.rm = TRUE) > 0)
metric_cols_lda_clean <- metric_cols_lda[keep]

lda_ref <- lda(Genotype ~ ., data = df_logo_ref[, c("Genotype", metric_cols_lda_clean)])
ref_scaling <- lda_ref$scaling[, lda_dim_names, drop = FALSE]  # reference loading vectors

# =========================================================
# LOGO loop — fit LDA leaving out each genotype in turn,
# collect sign-aligned loadings and left-out projections
# =========================================================
genotypes <- sort(unique(df_logo$Genotype))
n_geno    <- length(genotypes)

all_loadings       <- list()
pooled_projections <- list()

for (i in seq_along(genotypes)) {
  g_out <- genotypes[i]
  message("=== LOGO fold ", i, "/", n_geno, ": leaving out ", g_out, " ===")

  train_df <- df_logo %>% filter(Genotype != g_out)
  test_df  <- df_logo %>% filter(Genotype == g_out)

  # scale using training statistics only
  train_means <- sapply(train_df[metric_cols], mean, na.rm = TRUE)
  train_sds   <- sapply(train_df[metric_cols], sd,   na.rm = TRUE)
  train_sds[train_sds == 0] <- 1

  train_scaled <- train_df
  test_scaled  <- test_df
  train_scaled[metric_cols] <- scale(train_df[metric_cols], center = train_means, scale = train_sds)
  test_scaled[metric_cols]  <- scale(test_df[metric_cols],  center = train_means, scale = train_sds)

  # fit LDA without g_out
  lda_model <- lda(Genotype ~ ., data = train_scaled[, c("Genotype", metric_cols_lda_clean)])

  # sign-align each LD to the reference via dot product of loading vectors
  fold_scaling <- lda_model$scaling[, lda_dim_names, drop = FALSE]
  ld_signs <- sapply(lda_dim_names, function(d) {
    if (sum(ref_scaling[, d] * fold_scaling[, d]) < 0) -1L else 1L
  })

  # loadings (sign-aligned)
  ld_df <- data.frame(sweep(fold_scaling, 2, ld_signs, `*`))
  ld_df$Variable <- rownames(ld_df)
  ld_df$Left_out <- g_out
  all_loadings[[i]] <- ld_df

  # project left-out genotype onto fold LDA space; apply sign alignment
  lda_x_test <- data.frame(predict(lda_model, newdata = test_scaled)$x)[, lda_dim_names, drop = FALSE]
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

# =========================================================
# H2 of LD1-5 scores (pooled left-out projections)
# =========================================================
pooled_df <- bind_rows(pooled_projections)
h2_logo   <- calc_H2_fixedtime(pooled_df, lda_dim_names, day_col = "Day")
write.csv(h2_logo, result_path("logo_LDA_h2.csv"), row.names = FALSE)

p_h2_logo <- ggplot(h2_logo, aes(x = reorder(Trait, H2), y = H2)) +
  geom_col(fill = "salmon") +
  geom_text(aes(label = round(H2, 3)), hjust = -0.1, size = 3.5) +
  coord_flip() + ylim(0, max(h2_logo$H2) * 1.15) +
  theme_minimal() +
  labs(x = NULL, y = expression(H^2))
ggsave(result_path("logo_LDA_h2.png"), p_h2_logo, width = 6, height = 4, dpi = 300)

# =========================================================
# Mean loadings across LOGO folds
# =========================================================
loadings_all <- bind_rows(all_loadings)

mean_loadings <- loadings_all %>%
  group_by(Variable) %>%
  summarise(across(all_of(lda_dim_names), \(x) mean(x, na.rm = TRUE), .names = "mean_{.col}"),
            .groups = "drop")

write.csv(mean_loadings, result_path("logo_LDA_mean_loadings.csv"), row.names = FALSE)

# =========================================================
# Plots: top-20 mean LD1 / LD2 loadings, bar fill = raw trait H2
# (same style as make_loadplot() in 5_lda.R)
# =========================================================
raw_h2 <- read.csv(h2_path(sprintf("h2_timefixed_day_%s.csv", segmentation_method))) %>% select(Trait, H2)

load_df <- mean_loadings %>%
  left_join(raw_h2, by = c("Variable" = "Trait"))

make_loadplot <- function(ld) {
  col <- paste0("mean_", ld)
  top <- load_df %>% arrange(desc(abs(.data[[col]]))) %>% slice_head(n = 20)
  ggplot(top, aes(reorder(Variable, .data[[col]]), .data[[col]], fill = H2)) +
    geom_col() + coord_flip() +
    scale_fill_viridis(option = "plasma", name = expression(H^2)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 8)) +
    labs(x = NULL, y = sprintf("Mean loading across LOGO folds (%s)", ld))
}

ggsave(result_path("logo_LDA_loadings_LD1_h2color.png"), make_loadplot("LD1"), width = 6, height = 5.5, dpi = 300)
ggsave(result_path("logo_LDA_loadings_LD2_h2color.png"), make_loadplot("LD2"), width = 6, height = 5.5, dpi = 300)

cat("LOGO CV mean loadings and H2 complete. Outputs in", RESULT_DIR, "\n")
