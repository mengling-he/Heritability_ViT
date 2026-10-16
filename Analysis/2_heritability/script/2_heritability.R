# =============================================================================
# 2_heritability.R  — segmentation pipeline (U-Net / ViT), Step 2
# Broad-sense H2 for all side+top traits:
#   (a) per-date         : trait ~ (1|Genotype)                 [calc_H2_timeseries]
#   (b) date-as-covariate: trait ~ Day + (1|Genotype)+(1|G:Rep) [calc_H2_fixedtime]
# Produces the date-covariate H2 table + barplot, the per-date H2 table,
# and the timestamp-heritability trajectory plot (report Fig 4, Side/Top panels).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(lme4); library(ggplot2); library(stringr)
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

code_path <- function(...) here(base_dir, "Code", ...)
prep_path <- function(...) here(base_dir, "Analysis", "1_preprocessing", "result", ...)
result_path <- function(...) here(base_dir, "Analysis", "2_heritability", "result", ...)

RESULT_DIR <- result_path()
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

source(code_path("heritability_fun.R"))
source(code_path("plot_fun.R"))

# ---- helper: parse Genotype / Replicate from Plant.Info ----------------------
add_geno_rep <- function(df) {
  df %>% tidyr::extract(Plant.Info, into = c("Genotype", "Replicate"),
                        regex = "^\\d+_(.+)_(\\d+)$", remove = FALSE, convert = FALSE)
}

# =============================================================================
# Load daily-aligned data, drop slow growers, build modelling frame
# =============================================================================
df0 <- read.csv(prep_path(sprintf("side_top_features_%s_linearimputed_day.csv", segmentation_method)),
                stringsAsFactors = FALSE)
df0$date <- as.Date(df0$date)

slow <- strsplit(trimws(readLines(prep_path("slow_grower.txt"), warn = FALSE)), ",")[[1]]
df <- df0 %>% filter(!Plant.Info %in% slow) %>% add_geno_rep()
df$Day      <- as.numeric(df$date - min(df$date))
df$Date_cat <- factor(df$date)

cat(sprintf("Modelling frame: %d rows, %d plants, %d genotypes\n",
            nrow(df), n_distinct(df$Plant.Info), n_distinct(df$Genotype)))

# trait columns + drop near-constant traits (no variance -> singular / undefined H2)
trait_cols <- setdiff(names(df), c("Plant.Info", "date", "Genotype", "Replicate", "Day", "Date_cat"))
nzv <- sapply(df[trait_cols], function(x) {
  v <- var(x, na.rm = TRUE); is.finite(v) && v > .Machine$double.eps
})
dropped <- trait_cols[!nzv]
trait_cols <- trait_cols[nzv]
if (length(dropped)) cat("Dropped", length(dropped), "near-constant traits:",
                         paste(dropped, collapse = ", "), "\n")

# ---- keep only dates where ALL samples have values (no NA across modelled traits)
# Edge gaps from internal-only interpolation leave NA on early/late dates for
# plants whose first/last observation is not on the global min/max date. Drop
# those incomplete dates so every retained date is fully observed for every plant.
n_plants  <- n_distinct(df$Plant.Info)
date_ok <- df %>%
  group_by(date) %>%
  summarise(n_complete = sum(complete.cases(across(all_of(trait_cols)))),
            n_rows = n(), .groups = "drop") %>%
  filter(n_complete == n_rows & n_rows == n_plants)
keep_dates <- sort(date_ok$date)
dropped_dates <- setdiff(sort(unique(df$date)), keep_dates)

df <- df %>% filter(date %in% keep_dates)
df$Day      <- as.numeric(df$date - min(df$date))   # recompute after filtering
df$Date_cat <- factor(df$date)
cat(sprintf("Kept %d of %d dates (all samples present); dropped %d incomplete dates: %s\n",
            length(keep_dates), length(keep_dates) + length(dropped_dates),
            length(dropped_dates), paste(as.character(dropped_dates), collapse = ", ")))
cat(sprintf("Filtered modelling frame: %d rows, %d plants\n", nrow(df), n_distinct(df$Plant.Info)))

# =============================================================================
# (b) Date-as-covariate model:  trait ~ Day + (1|Genotype) + (1|Genotype:Replicate)
# =============================================================================
h2_fixed <- calc_H2_fixedtime(df, trait_cols, day_col = "Day")
write.csv(h2_fixed, result_path(sprintf("h2_timefixed_day_%s.csv", segmentation_method)), row.names = FALSE)

h2_fixed_sort <- h2_fixed %>% arrange(desc(H2))
p_fixed <- h2_fixed %>% filter(H2 > 0.2) %>%
  ggplot(aes(reorder(Trait, H2), H2)) +
  geom_col(fill = "#3B6EA5") +
  geom_text(aes(label = sprintf("%.2f", H2)), hjust = -0.1, size = 2.6) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(title = sprintf("Date-as-covariate H2 (%s, H2 > 0.2)", method_label), x = NULL, y = expression(H^2))
ggsave(result_path(sprintf("h2_timefixed_day_%s.png", segmentation_method)), p_fixed, width = 7, height = 8, dpi = 300)

# =============================================================================
# (a) Per-date model:  trait ~ (1|Genotype) at each observed date
# =============================================================================
h2_ts <- calc_H2_timeseries(df, trait_cols, time_col = "date")
h2_ts$Date <- as.Date(h2_ts$Time)
write.csv(h2_ts, result_path(sprintf("h2_timeseries_day_%s.csv", segmentation_method)), row.names = FALSE)

# heatmap of per-date H2 (traits with any H2 > 0.4)
h2_ts_high <- h2_ts %>% group_by(Trait) %>% filter(any(H2 > 0.4, na.rm = TRUE)) %>% ungroup()
plot_heatmap(h2_ts_high, columns_col = "Date", annot = FALSE,
             title = "Per-date heritability (H2 > 0.4 at some point)", show = FALSE,
             save_path = result_path(sprintf("h2_timeseries_heatmap_%s.png", segmentation_method)),
             figsize = c(14, 10))

# =============================================================================
# Report Fig 4: H2 trajectories over time, Side / Top panels
#   only traits with PEAK H2 > 0.45; strip view prefix for labels
# =============================================================================
h2_ts <- read.csv(result_path(sprintf("h2_timeseries_day_%s.csv", segmentation_method)),
                stringsAsFactors = FALSE)
peak_thresh <- 0.45
traj <- h2_ts %>%
  mutate(View  = ifelse(str_starts(Trait, "side_"), "Side View", "Top View"),
         Metric = str_remove(Trait, "^(side|top)_")) %>%
  group_by(Trait) %>% filter(max(H2, na.rm = TRUE) > peak_thresh) %>% ungroup()

cat(sprintf("Fig 4: %d traits with peak H2 > %.2f (%d side, %d top)\n",
            n_distinct(traj$Trait), peak_thresh,
            n_distinct(traj$Trait[traj$View == "Side View"]),
            n_distinct(traj$Trait[traj$View == "Top View"])))

p_traj <- ggplot(traj, aes(Date, H2, color = Metric, group = Metric)) +
  geom_hline(yintercept = peak_thresh, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  facet_wrap(~ View) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        legend.position = "right", legend.key.size = unit(0.4, "cm")) +
  labs(title = "Broad-sense heritability trajectories over time",
       x = "Date", y = expression(H^2), color = "Trait")
ggsave(result_path(sprintf("h2_trajectories_over_time_%s.png", segmentation_method)), p_traj,
       width = 13, height = 5.5, dpi = 300)

cat("\nStep 2 heritability complete. Outputs in", RESULT_DIR, "\n")


# --- Build one plot per view ---
plot_df <- h2_ts %>%
  mutate(View  = ifelse(str_starts(Trait, "side_"), "Side View", "Top View"),
         Metric = str_remove(Trait, "^(side|top)_")) %>%
  group_by(Trait) %>% filter(max(H2, na.rm = TRUE) > peak_thresh) %>% ungroup()


plot_df <- plot_df %>%
  group_by(View, Metric) %>%
  filter(any(H2 > 0.45)) %>%
  ungroup()
make_view_plot <- function(view_label) {
  ggplot(plot_df %>% filter(View == view_label),
         aes(x = Date, y = H2, color = Metric, group = Metric)) +
    ylim(0,0.7)+
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0.45, linetype = "dashed", color = "grey50") +
    labs(title = view_label,
         x = "Date", y = expression(H^2), color = "Trait") +
    theme_minimal(base_size = 11) +
    # Let patchwork handle the final positioning by keeping it standard here
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1,size = 6)
    )
}

p_side <- make_view_plot("Side View")
p_top  <- make_view_plot("Top View")
library(ggplot2)
library(patchwork)
p_by_view_series <- p_side + p_top + plot_layout(guides = "keep")
ggsave(result_path(sprintf("h2_trajectories_by_view_%s.png", segmentation_method)), p_by_view_series,
       width = 16, height = 5.5, dpi = 300)
print(p_by_view_series)
