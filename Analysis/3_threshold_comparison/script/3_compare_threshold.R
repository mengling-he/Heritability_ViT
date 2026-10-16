# =============================================================================
# 3_compare_threshold.R  — segmentation pipeline (U-Net / ViT), Step 3
# Compare date-as-covariate H2 from the segmentation pipeline against the thresholding
# method (raw_data/h2_timefixed_day_threshold.csv).
# Reproduces report Fig 3 (grouped bars, traits with H2 > 0.2 in at least one).
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(here) })

source(here("renv", "activate.R"))
segmentation_method <- "vit"  # Change to "unet" when needed
base_dir <- switch(
  segmentation_method,
  "unet" = "Heritability_Unet",
  "vit"  = "Heritability_ViT",
  stop("Unknown segmentation method")
)
method_label <- switch(segmentation_method, "unet" = "U-Net", "vit" = "ViT")

raw_path <- function(...) here(base_dir, "raw_data", ...)
h2_path <- function(...) here(base_dir, "Analysis", "2_heritability", "result", ...)
result_path <- function(...) here(base_dir, "Analysis", "3_threshold_comparison", "result", ...)

RESULT_DIR <- result_path()
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

seg   <- read.csv(h2_path(sprintf("h2_timefixed_day_%s.csv", segmentation_method))) %>% dplyr::select(Trait, H2) %>% rename(Seg = H2)
thres <- read.csv(h2_path("h2_timefixed_day_threshold.csv")) %>% dplyr::select(Trait, H2) %>% rename(Threshold = H2)

cmp <- inner_join(seg, thres, by = "Trait") %>%
  mutate(Diff = Seg - Threshold)
write.csv(cmp, result_path(sprintf("h2_%s_vs_threshold.csv", segmentation_method)), row.names = FALSE)

cat(sprintf("Traits compared: %d | %s higher: %d | Threshold higher: %d | mean diff: %.3f\n",
            nrow(cmp), method_label, sum(cmp$Diff > 0), sum(cmp$Diff < 0), mean(cmp$Diff)))

# ---- Fig 3: grouped bar, traits with H2 > 0.2 in at least one method ----------
keep <- cmp %>% filter(Seg > 0.2 | Threshold > 0.2)
long <- keep %>%
  pivot_longer(c(Seg, Threshold), names_to = "Dataset", values_to = "H2") %>%
  mutate(Dataset = recode(Dataset, Seg = method_label, Threshold = "Threshold"))
order_traits <- keep %>% arrange(desc(Seg)) %>% pull(Trait)
long$Trait <- factor(long$Trait, levels = order_traits)

p <- ggplot(long, aes(Trait, H2, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = setNames(c("#34B6C6", "#E8413A"), c("Threshold", method_label))) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        # axis.text.y = element_text(size = 12),
        # axis.title.x = element_text(size = 14),
        # axis.title.y = element_text(size = 14),
        # legend.text  = element_text(size = 12),
        legend.position = "right") +
  labs(#title = "Time-adjusted heritability: Threshold vs U-Net segmentation",
       #subtitle = "Traits with H2 > 0.2 in at least one pipeline",
       x = "Trait", y = expression(H^2), fill = "Dataset")
ggsave(result_path(sprintf("h2_compare_threshold_%s.png", segmentation_method)), p, width = 12, height = 5.5, dpi = 300)

# scatter (one point per trait) for the full comparison
p2 <- ggplot(cmp, aes(Threshold, Seg, label = Trait)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(color = "#E8413A", alpha = 0.7, size = 2) +
  theme_minimal(base_size = 12) +
  labs(title = sprintf("H2: %s vs Threshold (date-as-covariate)", method_label),
       x = expression(H^2 ~ "(Threshold)"), y = bquote(H^2 ~ "(" * .(method_label) * ")"))
ggsave(result_path(sprintf("h2_scatter_threshold_%s.png", segmentation_method)), p2, width = 6, height = 6, dpi = 300)

# # =============================================================================
# # Compare against the ViT date-as-covariate H2 from the Heritability_update run
# # =============================================================================
# VIT_PATH <- here("Heritability_update", "Analysis", "2_daily_imputed_analysis",
#                  "result", "h2_timefixed_day_vit.csv")
# 
# if (file.exists(VIT_PATH)) {
#   vit <- read.csv(VIT_PATH) %>% select(Trait, H2) %>% rename(ViT = H2)
# 
#   cmp_vit <- inner_join(unet, vit, by = "Trait") %>% mutate(Diff = U_Net - ViT)
#   write.csv(cmp_vit, result_path("h2_unet_vs_vit.csv"), row.names = FALSE)
#   cat(sprintf("U-Net vs ViT: %d traits | U-Net higher: %d | ViT higher: %d | mean diff: %.3f\n",
#               nrow(cmp_vit), sum(cmp_vit$Diff > 0), sum(cmp_vit$Diff < 0), mean(cmp_vit$Diff)))
# 
#   # pairwise U-Net vs ViT bar (traits H2 > 0.2 in at least one)
#   keep_v <- cmp_vit %>% filter(U_Net > 0.2 | ViT > 0.2)
#   long_v <- keep_v %>%
#     pivot_longer(c(U_Net, ViT), names_to = "Dataset", values_to = "H2") %>%
#     mutate(Dataset = recode(Dataset, U_Net = "U-Net"))
#   long_v$Trait <- factor(long_v$Trait, levels = keep_v %>% arrange(desc(U_Net)) %>% pull(Trait))
#   p_v <- ggplot(long_v, aes(Trait, H2, fill = Dataset)) +
#     geom_col(position = position_dodge(width = 0.8), width = 0.75) +
#     scale_fill_manual(values = c("ViT" = "#7A5BA6", "U-Net" = "#E8413A")) +
#     theme_minimal(base_size = 11) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) +
#     labs(title = "Time-adjusted heritability: ViT vs U-Net segmentation",
#          subtitle = "Traits with H2 > 0.2 in at least one pipeline",
#          x = "Trait", y = expression(H^2), fill = "Dataset")
#   ggsave(result_path("h2_compare_vit_unet.png"), p_v, width = 12, height = 5.5, dpi = 300)
# 
#   p_vs <- ggplot(cmp_vit, aes(ViT, U_Net, label = Trait)) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
#     geom_point(color = "#E8413A", alpha = 0.7, size = 2) +
#     theme_minimal(base_size = 12) +
#     labs(title = "H2: U-Net vs ViT (date-as-covariate)",
#          x = expression(H^2 ~ "(ViT)"), y = expression(H^2 ~ "(U-Net)"))
#   ggsave(result_path("h2_scatter_vit_unet.png"), p_vs, width = 6, height = 6, dpi = 300)
# 
#   # combined three-way comparison: Threshold vs ViT vs U-Net
#   three <- unet %>%
#     inner_join(thres, by = "Trait") %>%
#     inner_join(vit,   by = "Trait")
#   write.csv(three, result_path("h2_threshold_vit_unet.csv"), row.names = FALSE)
# 
#   keep3 <- three %>% filter(U_Net > 0.2 | Threshold > 0.2 | ViT > 0.2)
#   long3 <- keep3 %>%
#     pivot_longer(c(Threshold, ViT, U_Net), names_to = "Dataset", values_to = "H2") %>%
#     mutate(Dataset = factor(recode(Dataset, U_Net = "U-Net"),
#                             levels = c("Threshold", "ViT", "U-Net")))
#   long3$Trait <- factor(long3$Trait, levels = keep3 %>% arrange(desc(U_Net)) %>% pull(Trait))
#   p3 <- ggplot(long3, aes(Trait, H2, fill = Dataset)) +
#     geom_col(position = position_dodge(width = 0.85), width = 0.8) +
#     scale_fill_manual(values = c("Threshold" = "#34B6C6", "ViT" = "#7A5BA6", "U-Net" = "#E8413A")) +
#     theme_minimal(base_size = 11) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) +
#     labs(title = "Time-adjusted heritability: Threshold vs ViT vs U-Net",
#          subtitle = "Traits with H2 > 0.2 in at least one pipeline",
#          x = "Trait", y = expression(H^2), fill = "Segmentation")
#   ggsave(result_path("h2_compare_threshold_vit_unet.png"), p3, width = 13, height = 5.5, dpi = 300)
# } else {
#   cat("ViT reference not found at", VIT_PATH, "- skipping ViT comparison.\n")
# }
# 
# cat("Step 3 comparison complete. Outputs in", RESULT_DIR, "\n")
