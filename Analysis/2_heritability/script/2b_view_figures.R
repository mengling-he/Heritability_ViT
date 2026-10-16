# =============================================================================
# 2b_view_figures.R  — segmentation pipeline (U-Net / ViT), Step 2 (view-split figures)
# Two figures from the step-2 heritability results:
#   Fig A: date-as-covariate H2, side vs top grouped bars, sorted by per-trait mean
#   Fig B: per-date H2 trajectories, Side / Top panels (peak H2 > 0.45)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(patchwork)
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

result_path <- function(...) here(base_dir, "Analysis", "2_heritability", "result", ...)
RESULT_DIR <- result_path()
VIEW_COLORS <- c(side = "#66C2A5", top = "#FC8D62")

# =============================================================================
# Fig A: H2 (date-as-covariate) comparison, side vs top, sorted by mean value
# =============================================================================
h2 <- read.csv(result_path(sprintf("h2_timefixed_day_%s.csv", segmentation_method))) %>%
  mutate(View   = ifelse(str_starts(Trait, "side_"), "side", "top"),
         Metric = str_remove(Trait, "^(side|top)_"))

# keep only metrics with H2 > 0.1 in at least one view
keep_metrics <- h2 %>% group_by(Metric) %>%
  filter(max(H2, na.rm = TRUE) > 0.1) %>% pull(Metric) %>% unique()
h2 <- h2 %>% filter(Metric %in% keep_metrics)

metric_order <- h2 %>% group_by(Metric) %>%
  summarise(mean_h2 = mean(H2), .groups = "drop") %>%
  arrange(desc(mean_h2)) %>% pull(Metric)
h2$Metric <- factor(h2$Metric, levels = metric_order)

pA <- ggplot(h2, aes(Metric, H2, fill = View)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.7) +
  scale_fill_manual(values = VIEW_COLORS) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "plain")) +
  labs(title = expression("Heritability (" * H^2 * ") Comparison Sorted by Mean Value"),
       x = "Trait", y = expression(H^2), fill = "View")
ggsave(result_path(sprintf("h2_timefixed_day_%s_by_view.png", segmentation_method)), pA, width = 14, height = 7, dpi = 300)

# =============================================================================
# Fig B: per-date H2 trajectories, Side / Top panels, traits with peak H2 > 0.45
# =============================================================================
peak_thresh <- 0.45
ts <- read.csv(result_path(sprintf("h2_timeseries_day_%s.csv", segmentation_method))) %>%
  mutate(Date   = as.Date(Date),
         View   = ifelse(str_starts(Trait, "side_"), "Side View", "Top View"),
         Metric = str_remove(Trait, "^(side|top)_")) %>%
  group_by(Trait) %>% filter(max(H2, na.rm = TRUE) > peak_thresh) %>% ungroup()

make_panel <- function(view_label) {
  ggplot(ts %>% filter(View == view_label),
         aes(Date, H2, color = Metric, group = Metric)) +
    geom_hline(yintercept = peak_thresh, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
    scale_y_continuous(limits = c(0, 0.7)) +
    labs(title = view_label, x = "Date", y = expression(H^2), color = "Trait") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(hjust = 1, size = 14),
          # axis.title.x = element_text(size = 14),
          # axis.title.y = element_text(size = 14),
          # plot.title = element_text(size = 16),
          legend.text  = element_text(size = 14)
          )
}
pB <- make_panel("Side View") + make_panel("Top View") + plot_layout(guides = "keep")
ggsave(result_path(sprintf("h2_series_day_%s_by_view.png", segmentation_method)), pB, width = 20, height = 5, dpi = 300)

cat(sprintf("Fig A: %d traits (side+top bars). Fig B: %d side, %d top traits with peak H2 > %.2f.\n",
            n_distinct(h2$Metric),
            n_distinct(ts$Trait[ts$View == "Side View"]),
            n_distinct(ts$Trait[ts$View == "Top View"]), peak_thresh))
cat(sprintf("Saved h2_timefixed_day_%s_by_view.png and h2_series_day_%s_by_view.png\n",
            segmentation_method, segmentation_method))
