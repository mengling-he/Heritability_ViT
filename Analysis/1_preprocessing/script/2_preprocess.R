# =============================================================================
# 1_preprocess.R  — segmentation pipeline (U-Net / ViT), Step 1
# Combine side + top segmentation features and build two temporally-aligned
# datasets (slow growers: see 0_slow_grower.R):
#   (a) daily-aligned  : internal linear interpolation onto a daily grid
#   (b) date-group      : ~4-day bins, interpolation across bins
# Reproduces report Fig 2 (missing pattern).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(ggplot2)
  library(naniar)
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
method_label <- switch(segmentation_method, "unet" = "U-Net", "vit" = "ViT")
# raw feature file suffix and acquisition-date column differ between pipelines
raw_suffix <- switch(segmentation_method, "unet" = "", "vit" = "_vit")
date_col   <- switch(segmentation_method, "unet" = "Acquisition.Date", "vit" = "Analyse.Date")

RAW_DIR    <- here(base_dir, "raw_data")
RESULT_DIR <- here(base_dir, "Analysis", "1_preprocessing", "result")
DATA_DIR   <- here(base_dir, "Analysis", "1_preprocessing", "data")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Read raw side / top features
# =============================================================================
side0 <- read.csv(file.path(RAW_DIR, sprintf("5.25.2024_CBI.0010_York_Poplar_RGB.side_features%s.csv", raw_suffix)),
                  stringsAsFactors = FALSE, check.names = TRUE)
top0  <- read.csv(file.path(RAW_DIR, sprintf("5.25.2024_CBI.0010_York_Poplar_RGB.top_features%s.csv", raw_suffix)),
                  stringsAsFactors = FALSE, check.names = TRUE)

# trait columns = everything from "area" onwards (drop id / status / empty cols)
trait_start <- which(names(side0) == "area")
side_traits <- names(side0)[trait_start:ncol(side0)]
top_traits  <- names(top0)[which(names(top0) == "area"):ncol(top0)]

prep_view <- function(df, traits, prefix) {
  df %>%
    mutate(date = as.Date(.data[[date_col]])) %>%
    group_by(Plant.Info, date) %>%                       # >1 image/day -> middle point (median)
    summarise(across(all_of(traits), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
    rename_with(~ paste0(prefix, .x), all_of(traits))
}

side_df <- prep_view(side0, side_traits, "side_")
top_df  <- prep_view(top0,  top_traits,  "top_")

# combined, un-imputed (one row per plant-date that was actually observed)
combined <- full_join(side_df, top_df, by = c("Plant.Info", "date")) %>%
  arrange(Plant.Info, date)

trait_cols <- setdiff(names(combined), c("Plant.Info", "date"))
cat(sprintf("Combined: %d plant-date rows, %d plants, %d dates, %d traits\n",
            nrow(combined), n_distinct(combined$Plant.Info),
            n_distinct(combined$date), length(trait_cols)))

write.csv(combined, file.path(RESULT_DIR, sprintf("side_top_features_%s_final.csv", segmentation_method)), row.names = FALSE)

# =============================================================================
# 2. Slow growers: NOT computed here. The shared list (same plants for every
#    segmentation method) comes from 0_slow_grower.R -> result/slow_grower.txt
# =============================================================================

# =============================================================================
# 3. Missing pattern heatmap (report Fig 2) — on side area, daily grid
# =============================================================================
area_tab <- combined %>%
  dplyr::select(Plant.Info, date, side_area) %>%
  pivot_wider(names_from = date, values_from = side_area)
area_mat <- as.data.frame(area_tab)
rownames(area_mat) <- area_mat$Plant.Info
area_mat$Plant.Info <- NULL
area_mat <- area_mat[, order(as.Date(colnames(area_mat)))]
# sort rows by missingness pattern
area_mat <- area_mat[do.call(order, lapply(area_mat, function(x) !is.na(x))), ]

png(file.path(RESULT_DIR, "missing_pattern.png"), width = 22, height = 8, units = "in", res = 300)
print(
  vis_miss(area_mat, warn_large_data = FALSE, show_perc_col = TRUE) +
    labs(
      #title = "Missing pattern of the data", 
      y = "Plants") +
    theme(axis.text.x = element_text(size = 24, angle = 80, hjust = 1),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          legend.text  = element_text(size = 22))
)
dev.off()
# =============================================================================
# 4. Daily-aligned dataset: internal linear interpolation + end extrapolation
#    (rule = 2: leading/trailing gaps filled by carrying the first/last value)
# =============================================================================
# every calendar day between the first and last imaging date (not only observed dates)
daily_dates <- seq(min(combined$date), max(combined$date), by = "day")
full_grid <- expand_grid(Plant.Info = unique(combined$Plant.Info), date = daily_dates)

impute_daily <- full_grid %>%
  left_join(combined, by = c("Plant.Info", "date")) %>%
  arrange(Plant.Info, date) %>%
  group_by(Plant.Info) %>%
  mutate(across(all_of(trait_cols),
                ~ zoo::na.approx(.x, x = date, na.rm = FALSE, rule = 2))) %>%  # rule=2: extrapolate ends
  ungroup()

n_na_daily <- sum(is.na(impute_daily$side_area))
cat(sprintf("Daily-aligned: %d rows (%d plants x %d days); side_area NA remaining: %d\n",
            nrow(impute_daily), n_distinct(impute_daily$Plant.Info), length(daily_dates), n_na_daily))
write.csv(impute_daily, file.path(RESULT_DIR, sprintf("side_top_features_%s_linearimputed_day.csv", segmentation_method)),
          row.names = FALSE)

# =============================================================================
# 5. Date-group dataset: ~4-day bins, then interpolate across bins
# =============================================================================
bin_breaks <- as.Date(c("2024-05-26","2024-05-30","2024-06-04","2024-06-08","2024-06-12",
                        "2024-06-16","2024-06-20","2024-06-24","2024-06-28","2024-07-01"))
bin_labels <- c("May26-29","May30-Jun02","Jun04-07","Jun08-11","Jun12-15",
                "Jun16-19","Jun20-23","Jun24-27","Jun28-30")

combined_bin <- combined %>%
  mutate(Timepoint = cut(date, breaks = bin_breaks, labels = bin_labels,
                         right = FALSE, include.lowest = TRUE))

group_mean <- combined_bin %>%
  group_by(Plant.Info, Timepoint) %>%
  summarise(across(all_of(trait_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

bin_grid <- expand_grid(Plant.Info = unique(combined$Plant.Info),
                        Timepoint  = factor(bin_labels, levels = bin_labels))

impute_group <- bin_grid %>%
  left_join(group_mean, by = c("Plant.Info", "Timepoint")) %>%
  mutate(Day = as.integer(Timepoint)) %>%
  arrange(Plant.Info, Day) %>%
  group_by(Plant.Info) %>%
  mutate(across(all_of(trait_cols), ~ zoo::na.approx(.x, x = Day, na.rm = FALSE, rule = 2))) %>%
  ungroup()

cat(sprintf("Date-group: %d rows (%d plants x %d bins); side_area NA remaining: %d\n",
            nrow(impute_group), n_distinct(impute_group$Plant.Info), length(bin_labels),
            sum(is.na(impute_group$side_area))))
write.csv(impute_group, file.path(RESULT_DIR, sprintf("side_top_features_%s_linearimputed_bygroup.csv", segmentation_method)),
          row.names = FALSE)


# area-over-time (imputed daily) quick look
p_area <- ggplot(impute_daily, aes(date, side_area, group = Plant.Info)) +
  geom_line(alpha = 0.15, linewidth = 0.3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Side area over time (daily-aligned)", x = "Date", y = "Side area")
ggsave(file.path(RESULT_DIR, "side_area_over_time_imputed_day.png"), p_area, width = 9, height = 5, dpi = 300)

cat("\nStep 1 preprocessing complete. Outputs in", normalizePath(RESULT_DIR), "\n")
