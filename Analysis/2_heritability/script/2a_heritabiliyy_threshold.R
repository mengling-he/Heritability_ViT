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

# ----read data----



# =============================================================================
features_df0_top_side <- read.csv(prep_path("side_top_features_thres_linearimputed_day.csv"),
                stringsAsFactors = FALSE)
colnames(features_df0_top_side)

##### delete dead  #####
slow_growers <- strsplit(trimws(readLines(prep_path("slow_grower.txt"), warn = FALSE)), ",")[[1]]
features_df_top_side <- features_df0_top_side %>% filter(!Plant.Info %in% slow_growers) %>% add_geno_rep()
##### generate a Day column  #####
features_df_top_side$date <- as.Date(features_df_top_side$date)
features_df_top_side$Day <- as.numeric(features_df_top_side$date - min(features_df_top_side$date))
colnames(features_df_top_side)

#------- model: time as fixed effect  -------
trait_cols <- setdiff(names(features_df_top_side), c("Plant.Info", "date", "Genotype", "Replicate", "Day", "Date_cat"))
nzv <- sapply(features_df_top_side[trait_cols], function(x) {
  v <- var(x, na.rm = TRUE); is.finite(v) && v > .Machine$double.eps
})
dropped <- trait_cols[!nzv]
trait_cols <- trait_cols[nzv]
if (length(dropped)) cat("Dropped", length(dropped), "near-constant traits:",
                         paste(dropped, collapse = ", "), "\n")

h2_results <- calc_H2_fixedtime(features_df_top_side,trait_cols, day_col = "Day")

write.csv(h2_results, result_path("h2_timefixed_day_threshold.csv"), row.names = FALSE)

h2_results_sort <- h2_results[order(-h2_results$H2), ]
sum(h2_results$H2 > 0.2)
highH2list <- subset(h2_results, H2 > 0.2)$Trait

ggplot(h2_results, aes(x = H2, y = Trait)) +
  geom_col(fill = "skyblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.y = element_text(size=2)) +
  labs(title = "H2 for different metric based on daily imputed data", x = "", y = "H2")

ggplot(h2_results_sort[1:sum(h2_results$H2 > 0.2),], 
       aes(x = reorder(Trait, H2), y = H2)) +
  geom_col(fill = "skyblue") +
  coord_flip() +  # makes it horizontal and easier to read
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 12)
  ) +
  labs(
    title = "Top H2 (> 0.2) based on daily imputed data",
    x = "Trait",
    y = "H2"
  )




#------- model: across time  -------
h2_results_timeseries <- calc_H2_timeseries(features_df_top_side,trait_cols,time_col='date')
h2_results_timeseries$Date <- as.Date(h2_results_timeseries$Time)
write.csv(h2_results_timeseries, result_path("h2_timeseries_day_threshold.csv"), row.names = FALSE)

filtered_H2_trait_timeseries <- h2_results_timeseries %>%
  group_by(Trait) %>%
  filter(any(H2 > 0.4)) %>%
  ungroup()
filtered_H2_trait_timeseries %>%
  summarise(n_unique_traits = n_distinct(Trait))

plot_heatmap(filtered_H2_trait_timeseries,columns_col = "Date",annot=FALSE,title = 'Heritability (H²>0.4) over Time')

ggplot(h2_results_timeseries %>%
         group_by(Trait) %>%
         filter(any(H2 > 0.5)) %>%
         ungroup(), aes(x = Date, y = H2, color = Trait, group = Trait)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_vline(
    xintercept = as.Date('2024-06-19'),   # <-- adjust year as needed
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +annotate(
    "text",
    x = as.Date("2024-06-19"),
    y = 1.05,
    label = "June 19",
    angle = 0,
    vjust = -0.5,
    hjust = 1,
    size = 3.5
  ) +
  labs(
    title = "Heritability (H²>0.5) over Time",
    x = "Date",
    y = expression(H^2),
    color = "Trait"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

