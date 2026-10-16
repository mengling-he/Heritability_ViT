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
df0 <- read.csv(
  prep_path(sprintf("side_top_features_%s_linearimputed_day.csv",segmentation_method)),
  stringsAsFactors = FALSE
)
slow_growers <- strsplit(trimws(readLines(
    prep_path("slow_grower.txt"),
    warn = FALSE
  )),
  ","
)[[1]]
df <- df0 %>% filter(!Plant.Info %in% slow_growers) %>%
  tidyr::extract(Plant.Info, into = c("Genotype", "Replicate"),
                 regex = "^\\d+_(.+)_(\\d+)$", remove = FALSE, convert = FALSE)
df$date <- as.Date(df$date)
df$Day <- as.numeric(df$date - min(df$date))
# trait set: drop near-constant
traits <- grep("^(side|top)_", names(df), value = TRUE)
nzv <- sapply(df[traits], function(x) { v <- var(x, na.rm = TRUE); is.finite(v) && v > 1e-8 })
traits <- traits[nzv]

df$Genotype <- factor(df$Genotype)
df_scaled <- df
df_scaled[traits] <- scale(df[traits])
# guard: drop any column that became non-finite after scaling
cat(sprintf("LDA on %d traits, %d rows, %d genotypes\n", length(traits), nrow(df_scaled), nlevels(df_scaled$Genotype)))

# =============================================================================
# Fit LDA, extract LD scores + loadings
# =============================================================================
lda_model <- lda(as.formula(paste("Genotype ~", paste(traits, collapse = "+"))), data = df_scaled)
lda_x <- as.data.frame(predict(lda_model, df_scaled)$x)
LD_dims <- colnames(lda_x)[1:5]

LDA_df <- cbind(df_scaled[c("date", "Genotype", "Replicate", "Day")], lda_x[LD_dims])
LDA_df$Date_cat <- factor(LDA_df$date)

# H2 of LD scores: date-as-covariate model
h2_lda <- calc_H2_fixedtime(LDA_df, LD_dims, day_col = "Day")
write.csv(h2_lda, result_path(sprintf("h2_lda_%s.csv", segmentation_method)), row.names = FALSE)

p_h2lda <- ggplot(h2_lda, aes(reorder(Trait, H2), H2)) +
  geom_col(fill = "#3B6EA5") +
  geom_text(aes(label = sprintf("%.3f", H2)), hjust = -0.1, size = 3.5) +
  coord_flip() + ylim(0, max(h2_lda$H2) * 1.15) +
  theme_minimal(base_size = 12) +
  labs(title = "LDA heritability (date-as-covariate)", x = NULL, y = expression(H^2))
ggsave(result_path(sprintf("h2_lda_%s.png", segmentation_method)), p_h2lda, width = 6, height = 4, dpi = 300)
# 
# # =============================================================================
# # Leave-one-genotype-out CV heritability for LD1 (report headline ~0.66)
# # =============================================================================
# genos <- levels(df_scaled$Genotype)
# proj <- vector("list", length(genos))
# for (i in seq_along(genos)) {
#   g <- genos[i]
#   train <- df_scaled[df_scaled$Genotype != g, ]
#   test  <- df_scaled[df_scaled$Genotype == g, ]
#   m <- tryCatch(lda(as.formula(paste("Genotype ~", paste(traits, collapse = "+"))), data = droplevels(train)),
#                 error = function(e) NULL)
#   if (is.null(m)) next
#   px <- as.data.frame(predict(m, test)$x)[, 1:min(5, ncol(predict(m, test)$x)), drop = FALSE]
#   colnames(px) <- paste0("LD", seq_len(ncol(px)))
#   proj[[i]] <- cbind(test[c("date", "Genotype", "Replicate", "Day")], px)
# }
# logo_df <- bind_rows(proj)
# logo_dims <- grep("^LD", names(logo_df), value = TRUE)
# h2_logo <- calc_H2_fixedtime(logo_df, logo_dims, day_col = "Day")
# write.csv(h2_logo, result_path(sprintf("h2_lda_logo_cv_%s.csv", segmentation_method)), row.names = FALSE)
# cat("LOGO-CV LD heritability:\n"); print(h2_logo[, c("Trait", "H2")], row.names = FALSE)
# 
# # =============================================================================
# # Report Fig 6: top-20 LD1 / LD2 loadings, bar fill = raw trait H2
# # =============================================================================
# raw_h2 <- read.csv(h2_path(sprintf("h2_timefixed_day_%s.csv", segmentation_method))) %>% select(Trait, H2)
# load_df <- as.data.frame(lda_model$scaling) %>% rownames_to_column("Trait") %>%
#   left_join(raw_h2, by = "Trait")
# 
# make_loadplot <- function(ld) {
#   top <- load_df %>% arrange(desc(abs(.data[[ld]]))) %>% slice_head(n = 20)
#   ggplot(top, aes(reorder(Trait, .data[[ld]]), .data[[ld]], fill = H2)) +
#     geom_col() + coord_flip() +
#     scale_fill_viridis(option = "plasma", name = expression(H^2)) +
#     theme_minimal(base_size = 11) +
#     theme(axis.text.y = element_text(size = 8)) +
#     labs(title = sprintf("Top 20 Traits' Contributions to %s", ld),
#          x = NULL, y = sprintf("Contribution to %s", ld))
# }
# ggsave(result_path(sprintf("lda_loadings_LD1_%s.png", segmentation_method)), make_loadplot("LD1"), width = 6, height = 5.5, dpi = 300)
# ggsave(result_path(sprintf("lda_loadings_LD2_%s.png", segmentation_method)), make_loadplot("LD2"), width = 6, height = 5.5, dpi = 300)
# 
# cat("\nStep 5 LDA complete. Outputs in", RESULT_DIR, "\n")
