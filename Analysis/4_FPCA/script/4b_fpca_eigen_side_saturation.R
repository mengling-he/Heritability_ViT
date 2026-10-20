# =============================================================================
# 4b_fpca_eigen_side_area.R  — diagnose low FPCA heritability of side_area
# Runs FPCA for side_area (and side_solidity as a high-H2 contrast), saves the
# eigenfunctions (grid CSV + plot) and decomposes the FPC-score variance into
# genotype vs residual to explain the near-zero H2.
# =============================================================================

suppressPackageStartupMessages({
  library(fda); library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(lme4)
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
result_path <- function(...) here(base_dir, "Analysis", "4_FPCA", "result", ...)

RESULT_DIR <- result_path()
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

source(code_path("heritability_fun.R"))
source(fpca_code_path("tuning_nointer.R"))
source(fpca_code_path("pca_fun.R"))
source(fpca_code_path("pca_score.R"))
source(fpca_code_path("fpca_analysis_single.R"))

# ---- data: combined un-imputed, slow growers removed (same as step 4) --------
df0 <- read.csv(prep_path(sprintf("side_top_features_%s_final.csv", segmentation_method)), stringsAsFactors = FALSE)
df0$date <- as.Date(df0$date)
slow <- strsplit(trimws(readLines(prep_path("slow_grower.txt"), warn = FALSE)), ",")[[1]]
df_clean <- df0 %>% filter(!Plant.Info %in% slow) %>% arrange(Plant.Info, date)
df_clean$day_id <- as.numeric(difftime(df_clean$date, min(df_clean$date), units = "days"))

analyse_trait <- function(trait) {
  res <- run_fpca_trait(df_clean, trait = trait, id_col = "Plant.Info",
                        time_col = "day_id", date_col = "date",
                        order = 4, K_int = 12, J = 500, K_est = 2)

  # ---- eigenfunctions on the grid -> CSV ----
  tt  <- res$mean$tt
  phi <- res$fpca$phi_grid[, 1:2, drop = FALSE]
  eig <- data.frame(day = tt, mean_curve = res$mean$mu_hat,
                    phi1 = phi[, 1], phi2 = phi[, 2])
  write.csv(eig, result_path(sprintf("fpca_eigen_%s.csv", trait)), row.names = FALSE)

  # ---- eigenfunction plot (the function's own phi panel) ----
  ggsave(result_path(sprintf("fpca_eigen_%s.png", trait)),
         res$plots$phi, width = 11, height = 7, dpi = 300)
  
  #### plot 1 
  tlimit = range(df_clean$day_id)
  ylimit = range(df_clean[trait])
  
  t.all = sort(unique(df_clean$day_id))
  mu_df <- df_clean %>%
    group_by(day_id) %>%
    summarise(mu_emp = mean(.data[[trait]], na.rm = TRUE)) %>%
    arrange(day_id)
  mu.emp <- mu_df$mu_emp
  
  png(
    filename = result_path(sprintf("trajectory_plot_compare%s.png", trait)),
    width = 11, height = 8,units = "in", res = 300
  )
  par(mfrow = c(1, 2))
  # ---- Plot 1 ----
  plot(tlimit, ylimit, type = "n",
       xlab = "Day", ylab = trait,
       cex.axis = 0.8, cex.lab = 1)
  
  for (pid in unique(df_clean$Plant.Info)) {
    sub <- df_clean[df_clean$Plant.Info == pid, ]
    sub <- sub[order(sub$day_id), ]
    lines(sub$day_id, sub[[trait]], col = "gray70", lwd = 1)
  }
  
  lines(t.all, mu.emp, col = "green4", lwd = 2)
  
  legend("topleft",
         legend = c("Individual plants", "Empirical mean"),
         col = c("gray70", "green4"),
         lwd = c(1, 2),
         cex = 0.6, bty = "n")
  
  # ---- Plot 2 ----
  plot(tlimit, ylimit, type = "n",
       xlab = "Day", ylab = trait,
       cex.axis = 0.8, cex.lab = 1)
  
  lines(res$mean$tt,
        res$mean$mu_hat,
        col = "green4", lwd = 2)
  
  lines(t.all, mu.emp, col = "red", lwd = 2)
  
  legend("topleft",
         legend = c("Estimated", "Empirical"),
         col = c("green4", "red"),
         lwd = 2, cex = 0.6, bty = "n")
  
  dev.off()
  
  ###########################################################################################
  ######## Stage 3: Prediction of PC scores by BLUP ########
  ############################################################################################
  ### Principal component scores ###
  PC <- res$scores$PC
  phi.fun <- res$fpca$phi_grid
  phi.fun.2 <- phi.fun[,1:2]
  phi.pc <- PC %*% t(phi.fun.2)
  dim(phi.pc)
  mu_vec <- res$mean$mu_hat
  y_hat  <- t(t(phi.pc) + mu_vec)
  
  tt <- res$mean$tt
  png(
    filename = result_path(sprintf("fitted_sample_%s.png", trait)),
    width = 11, height = 8,units = "in", res = 300
  )
  par(mfrow = c(2, 2))
  for (id in c(2,8,32,64)){
    id_ex <- rownames(PC)[id]
    sub_ex <- df_clean[df_clean$Plant.Info == id_ex, ]
    
    plot(sub_ex$day_id, sub_ex[[trait]],
         col="gray70", lwd=1,
         xlim = range(tt),xlab = "Day", ylab = trait,main = id_ex)
    lines(tt, y_hat[id,], col="blue", lwd=2)
  }
  dev.off()
  # pc_sub <- PC[c(1,399,229,363),]
  # colnames(pc_sub) <- c("fPC1_score", "fPC2_score", "fPC3_score")
  # row.names(pc_sub) <- names(Tlist)[c(1,399,229,363)]
  
  # ---- variance decomposition of FPC scores (genotype vs residual) ----
  pc <- res$scores$PC_df %>%
    tidyr::extract(Plant.Info, into = c("Genotype", "Replicate"),
                   regex = "^\\d+_(.+)_(\\d+)$", remove = FALSE, convert = FALSE)
  dec <- lapply(c("fpca1", "fpca2"), function(s) {
    m  <- lmer(as.formula(paste(s, "~ (1|Genotype)")), data = pc)
    vc <- as.data.frame(VarCorr(m))
    vg <- vc$vcov[vc$grp == "Genotype"]; ve <- vc$vcov[vc$grp == "Residual"]
    data.frame(score = s, var_genotype = vg, var_residual = ve,
               H2 = round(vg / (vg + ve), 4), score_sd = sd(pc[[s]]))
  })
  dec <- do.call(rbind, dec)

  list(trait = trait,
       pve = round(res$fpca$pve[1:3], 4),
       eigenvalues = round(res$fpca$vpos[1:3], 3),
       sigma_e2 = res$fpca$sigma_e2,
       decomp = dec)
}

diag_traits <- commandArgs(trailingOnly = TRUE)
if (length(diag_traits) == 0) diag_traits <- c("top_saturation_mean")
for (tr in diag_traits) {
  cat("\n==================== ", tr, " ====================\n")
  out <- analyse_trait(tr)
  cat("PVE (PC1-3):", paste(out$pve, collapse = ", "), "\n")
  cat("Eigenvalues:", paste(out$eigenvalues, collapse = ", "), "\n")
  cat("sigma_e2   :", format(out$sigma_e2, scientific = TRUE, digits = 3), "\n")
  cat("FPC-score variance decomposition (genotype vs residual):\n")
  print(out$decomp, row.names = FALSE)
}
cat("\nEigenfunction CSVs and plots saved to", RESULT_DIR, "\n")








