# =============================================================================
# 4_fpca.R  — segmentation pipeline (U-Net / ViT), Step 4
# Functional PCA on smoothed trait trajectories (un-imputed combined data).
# For each trait: estimate FPC1/FPC2 scores per plant, then broad-sense H2 of
# those scores (trait ~ (1|Genotype), calc_H2_1dim).
# Reproduces report Fig 5 (grouped bars FPC1 / FPC2 / LMM_w_time).
# =============================================================================

suppressPackageStartupMessages({
  library(fda); library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(here)
})

source(here("renv", "activate.R"))
segmentation_method <- "vit"  # Change to "unet" when needed
base_dir <- switch(
  segmentation_method,
  "unet" = "Heritability_Unet",
  "vit"  = "Heritability_ViT",
  stop("Unknown segmentation method")
)
code_path <- function(...) here(base_dir, "Code", ...)
fpca_code_path <- function(...) here(base_dir, "Code", "FPCA_code", ...)
prep_path <- function(...) here(base_dir, "Analysis", "1_preprocessing", "result", ...)
h2_path <- function(...) here(base_dir, "Analysis", "2_heritability", "result", ...)
result_path <- function(...) here(base_dir, "Analysis", "4_FPCA", "result2", ...)
RESULT_DIR <- result_path()
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
source(code_path("heritability_fun.R"))
source(fpca_code_path("fpca_analysis_single.R"))
# =============================================================================
# Data: combined un-imputed trajectories, slow growers removed
# =============================================================================
df0 <- read.csv(prep_path(sprintf("side_top_features_%s_final.csv", segmentation_method)), stringsAsFactors = FALSE)
df0$date <- as.Date(df0$date)
slow <- strsplit(trimws(readLines(prep_path("slow_grower.txt"), warn = FALSE)), ",")[[1]]
df_clean <- df0 %>% filter(!Plant.Info %in% slow) %>% arrange(Plant.Info, date)
df_clean$day_id <- as.numeric(difftime(df_clean$date, min(df_clean$date), units = "days"))

# traits: side_/top_, drop near-constant (FPCA needs variation)
traits <- grep("^(side|top)_", names(df_clean), value = TRUE)
nzv <- sapply(df_clean[traits], function(x) { v <- var(x, na.rm = TRUE); is.finite(v) && v > 1e-8 })
traits <- traits[nzv]

# Standardize each trait (z-score) before FPCA. Large-magnitude traits (e.g.
# side_area ~1e5-1e6) otherwise produce FPC scores ~1e11 whose variance-component
# fit is singular, collapsing genotype variance to 0 and giving a spurious H2~0.
# H2 is scale-invariant in theory; standardizing fixes the numerical conditioning.
df_clean <- df_clean %>% mutate(across(all_of(traits), ~ as.numeric(scale(.x))))
cat("FPCA on", length(traits), "traits (standardized)\n")

# =============================================================================
# Run FPCA per trait, collect FPC1/FPC2 scores and their H2
# =============================================================================
fpc1_list <- list(); fpc2_list <- list()
h2_rows   <- list()

for (y in traits) {
  res <- tryCatch(
    run_fpca_trait(df_clean, trait = y, id_col = "Plant.Info",
                   time_col = "day_id", date_col = "date",
                   order = 4, K_int = 12, J = 500, K_est = 2),
    error = function(e) { message("FPCA failed for ", y, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) next

  pc <- res$scores$PC_df %>%
    tidyr::extract(Plant.Info, into = c("Genotype", "Replicate"),
                   regex = "^\\d+_(.+)_(\\d+)$", remove = FALSE, convert = FALSE)

  avail <- intersect(c("fpca1", "fpca2"), names(pc))
  H2 <- tryCatch(calc_H2_1dim(pc, avail), error = function(e) NULL)
  if (is.null(H2)) next

  h2_rows[[y]] <- tibble(Trait = y,
                         FPC1 = H2$H2[match("fpca1", H2$Trait)],
                         FPC2 = H2$H2[match("fpca2", H2$Trait)])

  if ("fpca1" %in% names(pc))
    fpc1_list[[y]] <- pc %>% dplyr::select(Plant.Info, !!y := fpca1)
  if ("fpca2" %in% names(pc))
    fpc2_list[[y]] <- pc %>% dplyr::select(Plant.Info, !!y := fpca2)
}

fpc_h2 <- bind_rows(h2_rows)
write.csv(fpc_h2, result_path("fpc_h2_all.csv"), row.names = FALSE)

join_scores <- function(lst) Reduce(function(a, b) full_join(a, b, by = "Plant.Info"), lst)
if (length(fpc1_list)) write.csv(join_scores(fpc1_list), result_path("fpc1_score_all.csv"), row.names = FALSE)
if (length(fpc2_list)) write.csv(join_scores(fpc2_list), result_path("fpc2_score_all.csv"), row.names = FALSE)

# =============================================================================
# Report Fig 5: FPC1 / FPC2 / LMM_w_time for union of top-10 traits of each
# =============================================================================
lmm <- read.csv(h2_path(sprintf("h2_timefixed_day_%s.csv",segmentation_method)),
                stringsAsFactors = FALSE) %>%dplyr::select(Trait, LMM_w_time = H2)
fpc_h2 <-  read.csv(result_path("fpc_h2_all.csv"))
merged <- fpc_h2 %>% left_join(lmm, by = "Trait")

top10 <- function(col) merged %>% arrange(desc(.data[[col]])) %>% slice_head(n = 10) %>% pull(Trait)
sel <- unique(c(top10("FPC1"), top10("FPC2"), top10("LMM_w_time")))
sel <- sel[!is.na(sel)]

plot_df <- merged %>% filter(Trait %in% sel) %>%
  pivot_longer(c(FPC1, FPC2, LMM_w_time), names_to = "Metric", values_to = "H2") %>%
  mutate(Metric = factor(Metric, levels = c("FPC1", "FPC2", "LMM_w_time")))
order_traits <- merged %>% filter(Trait %in% sel) %>% arrange(desc(LMM_w_time)) %>% pull(Trait)
plot_df$Trait <- factor(plot_df$Trait, levels = order_traits)

p <- ggplot(plot_df, aes(Trait, H2, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("FPC1" = "#2E8B57", "FPC2" = "#7FBF7F", "LMM_w_time" = "#C0392B")) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14)) +
  labs(#title = "Heritability estimates from FPCA and LMM with time effect",
       x = "Trait", y = "Heritability", fill = "Metric")
ggsave(result_path(sprintf("h2_FPCA_vs_LMM_%s.png",segmentation_method)), p, width = 11, height = 5.5, dpi = 300)

cat(sprintf("FPCA H2 computed for %d traits; Fig 5 uses %d union traits.\n", nrow(fpc_h2), length(sel)))
cat("Step 4 FPCA complete. Outputs in", RESULT_DIR, "\n")
