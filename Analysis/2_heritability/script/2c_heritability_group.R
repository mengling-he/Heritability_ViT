# =============================================================================
# 2c_heritability_group.R  — segmentation pipeline (U-Net / ViT), Step 2
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
# =============================================================================
# Load daily-aligned data, drop slow growers, build modelling frame
# =============================================================================
df0 <- read.csv(
  prep_path(sprintf("side_top_features_%s_linearimputed_bygroup.csv",segmentation_method)),
                stringsAsFactors = FALSE)
df0$Day <- factor(df0$Day)

slow <- strsplit(trimws(readLines(
  prep_path("slow_grower.txt"),
  warn = FALSE
)),
","
)[[1]]
df <- df0 %>% filter(!Plant.Info %in% slow) %>% 
  tidyr::extract(Plant.Info, into = c("Genotype", "Replicate"),regex = "^\\d+_(.+)_(\\d+)$", 
                 remove = FALSE, convert = FALSE)
cat(sprintf("Modelling frame: %d rows, %d plants, %d genotypes\n",
            nrow(df), n_distinct(df$Plant.Info), n_distinct(df$Genotype)))

# trait columns + drop near-constant traits (no variance -> singular / undefined H2)
trait_cols <- setdiff(names(df), c("Plant.Info", "Timepoint", "Genotype", "Replicate", "Day"))
nzv <- sapply(df[trait_cols], function(x) {
  v <- var(x, na.rm = TRUE); is.finite(v) && v > .Machine$double.eps
})
dropped <- trait_cols[!nzv]
trait_cols <- trait_cols[nzv]
if (length(dropped)) cat("Dropped", length(dropped), "near-constant traits:",
                         paste(dropped, collapse = ", "), "\n")

# ---- keep only dates where ALL samples have values (no NA across modelled traits)
n_plants  <- n_distinct(df$Plant.Info)
date_ok <- df %>%
  group_by(Timepoint) %>%
  summarise(n_complete = sum(complete.cases(across(all_of(trait_cols)))),
            n_rows = n(), .groups = "drop") %>%
  filter(n_complete == n_rows & n_rows == n_plants)
keep_dates <- sort(date_ok$Timepoint)
dropped_dates <- setdiff(sort(unique(df$Timepoint)), keep_dates)

df <- df %>% filter(Timepoint %in% keep_dates)
cat(sprintf("Kept %d of %d dates (all samples present); dropped %d incomplete dates: %s\n",
            length(keep_dates), length(keep_dates) + length(dropped_dates),
            length(dropped_dates), paste(as.character(dropped_dates), collapse = ", ")))
cat(sprintf("Filtered modelling frame: %d rows, %d plants\n", nrow(df), n_distinct(df$Plant.Info)))

# =============================================================================
# (a) Date-as-covariate model:  trait ~ Day + (1|Genotype) + (1|Genotype:Replicate)
# =============================================================================
h2_fixed <- calc_H2_fixedtime(df, trait_cols, day_col = "Day")
write.csv(h2_fixed, result_path(sprintf("h2_timefixed_bygroup_%s.csv",segmentation_method)), row.names = FALSE)

h2_fixed_sort <- h2_fixed %>% arrange(desc(H2))
p_fixed <- h2_fixed %>% filter(H2 > 0.2) %>%
  ggplot(aes(reorder(Trait, H2), H2)) +
  geom_col(fill = "#3B6EA5") +
  geom_text(aes(label = sprintf("%.2f", H2)), hjust = -0.1, size = 2.6) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(title = sprintf("Date-as-covariate H2 (%s, H2 > 0.2)", method_label), x = NULL, y = expression(H^2))
ggsave(result_path(sprintf("h2_timefixed_bygroup_%s.png",segmentation_method)), p_fixed, width = 7, height = 8, dpi = 300)


# =============================================================================
# (b) compare with daily
# =============================================================================
h2_fixed_day  <- read.csv(result_path(sprintf("h2_timefixed_day_%s.csv",segmentation_method)))%>% 
  dplyr::select(Trait, H2) %>% rename(byDay = H2)
h2_fixed_group  <- read.csv(result_path(sprintf("h2_timefixed_bygroup_%s.csv",segmentation_method)))%>% 
  dplyr::select(Trait, H2) %>% rename(byGroup = H2)
cmp <- inner_join(h2_fixed_day, h2_fixed_group, by = "Trait") %>%
  mutate(Diff = byDay - byGroup)
write.csv(cmp, result_path(sprintf("h2_%s_imputation.csv", segmentation_method)), row.names = FALSE)
keep <- cmp %>% filter(byDay > 0.2 | byGroup > 0.2)
long <- keep %>%
  pivot_longer(c(byDay, byGroup), names_to = "Dataset", values_to = "H2") %>%
  mutate(Dataset = recode(Dataset, byDay = "byDay", byGroup = "byGroup"))
order_traits <- keep %>% arrange(desc(byDay)) %>% pull(Trait)
long$Trait <- factor(long$Trait, levels = order_traits)

p <- ggplot(long, aes(Trait, H2, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("byGroup" = "#34B6C6", "byDay" = "#E8413A")) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        # axis.text.y = element_text(size = 12),
        # axis.title.x = element_text(size = 14),
        # axis.title.y = element_text(size = 14),
        # legend.text  = element_text(size = 12),
        legend.position = "right") +
  labs(#title = "Time-adjusted heritability: byGroup vs byDay Imputation",
    #subtitle = "Traits with H2 > 0.2 in at least one pipeline",
    x = "Trait", y = expression(H^2), fill = "Dataset")
ggsave(result_path("h2_compare_imputation.png"), p, width = 12, height = 5.5, dpi = 300)
