# =============================================================================
# 0_slow_grower.R  — Step 0 (run before 1_preprocess.R / 1a_preprocess_thresholding.R)
# Shared slow-grower (non-viable plant) list used by EVERY analysis, so all
# segmentation methods (thresholding / ViT / U-Net) are compared on the same plants.
# Rule: side-area growth < 5% between a plant's first and last imaging date,
# computed on the THRESHOLDING features (this reproduces the 13-plant list used
# in Heritability_update).
# Output: 1_preprocessing/result/slow_grower.txt
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(here)
})

source(here("renv", "activate.R"))
segmentation_method <- "vit"  # Change to "unet" when needed (only sets the output folder)
base_dir <- switch(
  segmentation_method,
  "unet" = "Heritability_Unet",
  "vit"  = "Heritability_ViT",
  stop("Unknown segmentation method")
)

RAW_DIR    <- here(base_dir, "raw_data")
RESULT_DIR <- here(base_dir, "Analysis", "1_preprocessing", "result")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Thresholding side / top features -> one row per plant-date (median)
# =============================================================================
side0 <- read.csv(file.path(RAW_DIR, "5.25.2024_CBI.0010_York_Poplar_RGB.side_features_thresholding.csv"),
                  stringsAsFactors = FALSE, check.names = TRUE)
top0  <- read.csv(file.path(RAW_DIR, "5.25.2024_CBI.0010_York_Poplar_RGB.top_features_thresholding.csv"),
                  stringsAsFactors = FALSE, check.names = TRUE)

prep_view <- function(df, prefix) {
  df %>%
    mutate(date = as.Date(Analyse.Date)) %>%
    group_by(Plant.Info, date) %>%                       # >1 image/day -> middle point (median)
    summarise(area = median(area, na.rm = TRUE), .groups = "drop") %>%
    rename_with(~ paste0(prefix, .x), "area")
}

# same side+top date union as step 1, so first/last dates match
combined <- full_join(prep_view(side0, "side_"), prep_view(top0, "top_"),
                      by = c("Plant.Info", "date")) %>%
  arrange(Plant.Info, date)

# =============================================================================
# 2. Slow growers: side-area growth < 5% over the season
# =============================================================================
slow_calc <- combined %>%
  group_by(Plant.Info) %>%
  arrange(date, .by_group = TRUE) %>%
  summarise(
    first_area = side_area[which.min(date)],
    last_area  = side_area[which.max(date)],
    growth_pct = (last_area - first_area) / first_area * 100,
    .groups = "drop"
  ) %>%
  filter(growth_pct < 5)

slow <- sort(unique(slow_calc$Plant.Info))
writeLines(paste(slow, collapse = ","), file.path(RESULT_DIR, "slow_grower.txt"))
cat(sprintf("Slow growers (thresholding side area, growth < 5%%): %d plants\n", length(slow)))
print(as.data.frame(slow_calc), row.names = FALSE)

# check against the reference list from Heritability_update
ref_path <- here("Heritability_update", "Analysis", "1_preprocessing", "data", "slow_grower.txt")
if (file.exists(ref_path)) {
  slow_ref <- sort(strsplit(trimws(readLines(ref_path, warn = FALSE)), ",")[[1]])
  cat("Identical to Heritability_update slow_grower.txt:", setequal(slow, slow_ref), "\n")
}

# trajectories of the slow growers (sanity check plot)
p_slow <- ggplot(combined %>% filter(Plant.Info %in% slow),
                 aes(date, side_area, group = Plant.Info, color = Plant.Info)) +
  geom_line() + geom_point(size = 1) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Slow growers (thresholding side area over time)", x = "Date", y = "Side area")
ggsave(file.path(RESULT_DIR, "slow_grower_trajectories.png"), p_slow, width = 9, height = 5, dpi = 300)

cat("Saved", file.path(RESULT_DIR, "slow_grower.txt"), "\n")
