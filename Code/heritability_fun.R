# =============================================================================
# heritability.R
# Broad-sense heritability (H²) estimation and bootstrapping
#
# Functions:
#   calc_H2_1dim            — single H² across all data (no time)
#   calc_H2_fixedtime       — H² with Day as fixed effect
#   calc_H2_timeseries      — H² per time point
#   bootstrap_H2_1dim       — bootstrap CI for calc_H2_1dim
#   bootstrap_H2_fixedtime  — bootstrap CI for calc_H2_fixedtime
#   bootstrap_H2_timeseries — bootstrap CI per time point
# =============================================================================

library(lme4)
library(dplyr)
library(tidyr)
library(purrr)

# =============================================================================
# Internal helpers (not exported)
# =============================================================================

#' Fit an lmer model safely; returns NULL on failure
.fit_lmer <- function(trait, formula_str, data) {
  tryCatch(
    lmer(as.formula(paste(trait, formula_str)), data = data),
    error = function(e) {
      warning("Model failed for trait '", trait, "': ", conditionMessage(e))
      NULL
    }
  )
}

#' Extract variance components as a 1-row tibble; fills missing terms with 0
.extract_varcor <- function(model,
                            components = c("Genotype", "Genotype:Replicate", "Residual")) {
  vc <- as.data.frame(VarCorr(model)) %>%
    dplyr::select(grp, vcov) %>%
    tidyr::pivot_wider(names_from = grp, values_from = vcov)
  
  missing <- setdiff(components, colnames(vc))
  vc[missing] <- 0
  vc
}

#' Compute broad-sense H² (rounded to 4 dp)
.calc_h2 <- function(sigma_G2, sigma_E2) {
  round(sigma_G2 / (sigma_G2 + sigma_E2), 4)
}

#' Run one bootstrap replicate and return a single H² value
#'
#' @param formula_str  RHS string e.g. "~ Day + (1|Genotype) + (1|Genotype:Replicate)"
#' @param vc_names     variance component names; first = G, rest sum to E
.boot_one <- function(data, trait, formula_str,
                      vc_names = c("Genotype", "Genotype:Replicate", "Residual")) {
  # Stratified bootstrap: resample within each genotype
  boot_data <- do.call(rbind, lapply(split(data, data$Genotype), function(g) {
    g[sample(nrow(g), replace = TRUE), ]
  }))
  model     <- .fit_lmer(trait, formula_str, boot_data)
  if (is.null(model)) return(NA_real_)
  
  vc       <- .extract_varcor(model, vc_names)
  sigma_G2 <- vc[[vc_names[1]]]
  sigma_E2 <- rowSums(vc[vc_names[-1]])   # all non-G components sum to E
  
  .calc_h2(sigma_G2, sigma_E2)
}

#' Summarise a vector of bootstrap H² replicates into mean / SE / 95 % CI
.summarise_boot <- function(h2_boot) {
  ci <- quantile(h2_boot, probs = c(0.025, 0.975), na.rm = TRUE)
  list(
    h2_boot  = h2_boot,
    mean     = mean(h2_boot, na.rm = TRUE),
    se       = sd(h2_boot,   na.rm = TRUE),
    ci_lower = ci[[1]],
    ci_upper = ci[[2]]
  )
}

# =============================================================================
# H² estimation
# =============================================================================

#' Broad-sense H² across all data (no time structure)
#'
#' Model: trait ~ (1|Genotype)
#'
#' @param data    data.frame
#' @param traits  character vector of trait column names
#' @return tibble with columns: Trait, Var_Genotype, Var_Residual, H2
calc_H2_1dim <- function(data, traits) {
  
  purrr::map_dfr(traits, function(trait) {
    model <- .fit_lmer(trait, "~ (1|Genotype)", data)
    if (is.null(model)) return(NULL)
    
    vc <- .extract_varcor(model, c("Genotype", "Residual"))
    
    tibble(
      Trait        = trait,
      Var_Genotype = vc$Genotype,
      Var_Residual = vc$Residual,
      H2           = .calc_h2(vc$Genotype, vc$Residual)
    )
  })
}

#' Broad-sense H² with Day as a fixed effect (single model across time)
#'
#' Model: trait ~ Day + (1|Genotype) + (1|Genotype:Replicate)
#'
#' @param data     data.frame
#' @param traits   character vector of trait column names
#' @param day_col  name of the time/day column (default "Day")
#' @return tibble with columns: Trait, Var_Genotype, Var_Replicate, Var_Residual, H2
calc_H2_fixedtime <- function(data, traits, day_col = "Day") {
  
  formula_str <- paste("~", day_col, "+ (1|Genotype) + (1|Genotype:Replicate)")
  
  purrr::map_dfr(traits, function(trait) {
    message("Processing: ", trait)
    model <- .fit_lmer(trait, formula_str, data)
    if (is.null(model)) return(NULL)
    
    vc       <- .extract_varcor(model, c("Genotype", "Genotype:Replicate", "Residual"))
    sigma_G2 <- vc$Genotype
    sigma_E2 <- vc$`Genotype:Replicate` + vc$Residual
    
    tibble(
      Trait         = trait,
      Var_Genotype  = sigma_G2,
      Var_Replicate = vc$`Genotype:Replicate`,
      Var_Residual  = vc$Residual,
      H2            = .calc_h2(sigma_G2, sigma_E2)
    )
  })
}

#' Broad-sense H² estimated separately for each time point
#'
#' Model per time point: trait ~ (1|Genotype)
#'
#' @param data         data.frame
#' @param traits       character vector of trait column names
#' @param time_col     name of the time/day column (default "Day")
#' @param time_groups  optional vector of time points; defaults to all unique values
#' @return tibble with columns: Time, Trait, Var_Genotype, Var_Residual, H2
calc_H2_timeseries <- function(data, traits, time_col = "Day", time_groups = NULL) {
  
  time_groups <- time_groups %||% sort(unique(data[[time_col]]))
  
  purrr::map_dfr(time_groups, function(d) {
    message("Processing Day: ", d)
    df_day <- dplyr::filter(data, .data[[time_col]] == d)
    
    purrr::map_dfr(traits, function(trait) {
      cat("  Trait:", trait, "\n")
      model <- .fit_lmer(trait, "~ (1|Genotype)", df_day)
      if (is.null(model)) return(NULL)
      
      vc <- .extract_varcor(model, c("Genotype", "Residual"))
      
      tibble(
        Time         = d,
        Trait        = trait,
        Var_Genotype = vc$Genotype,
        Var_Residual = vc$Residual,
        H2           = .calc_h2(vc$Genotype, vc$Residual)
      )
    })
  })
}

# =============================================================================
# Bootstrap H² (SE and 95 % CI via percentile method)
# =============================================================================

#' Bootstrap H² — no time (matches calc_H2_1dim)
#'
#' @param data   data.frame
#' @param trait  single trait name (character)
#' @param B      number of bootstrap replicates (default 1000)
#' @param seed   random seed for reproducibility (default 123)
#' @return list: h2_boot, mean, se, ci_lower, ci_upper
bootstrap_H2_1dim <- function(data, traits, B = 100, seed = 123) {
  set.seed(seed)
  formula_str <- "~ (1|Genotype) "
  vc_names    <- c("Genotype", "Residual")
  
  # Run bootstrap for every trait
  boot_list <- purrr::map(traits, function(trait) {
    message("Bootstrapping: ", trait)
    replicate(B, .boot_one(data, trait, formula_str, vc_names))
  }) %>% purrr::set_names(traits)
  
  # Summarise into a tidy dataframe
  summary_df <- purrr::map_dfr(traits, function(trait) {
    s <- .summarise_boot(boot_list[[trait]])
    tibble(
      Trait    = trait,
      H2_mean  = s$mean,
      H2_se    = s$se,
      CI_lower = s$ci_lower,
      CI_upper = s$ci_upper
    )
  })
  
  list(
    boot_distributions = boot_list,   # named list: trait -> numeric vector length B
    summary            = summary_df   # one row per trait
  )
}

#' Bootstrap H² — Day as fixed effect (matches calc_H2_fixedtime)
#'
#' @param data     data.frame
#' @param trait    single trait name (character)
#' @param day_col  name of the time/day column (default "Day")
#' @param B        number of bootstrap replicates (default 1000)
#' @param seed     random seed for reproducibility (default 123)
#' @return list: h2_boot, mean, se, ci_lower, ci_upper
bootstrap_H2_fixedtime <- function(data, traits, day_col = "Day", B = 100, seed = 123) {
  set.seed(seed)
  formula_str <- paste("~", day_col, "+ (1|Genotype) + (1|Genotype:Replicate)")
  vc_names    <- c("Genotype", "Genotype:Replicate", "Residual")
  
  # Run bootstrap for every trait
  boot_list <- purrr::map(traits, function(trait) {
    message("Bootstrapping: ", trait)
    replicate(B, .boot_one(data, trait, formula_str, vc_names))
  }) %>% purrr::set_names(traits)
  
  # Summarise into a tidy dataframe
  summary_df <- purrr::map_dfr(traits, function(trait) {
    s <- .summarise_boot(boot_list[[trait]])
    tibble(
      Trait    = trait,
      H2_mean  = s$mean,
      H2_se    = s$se,
      CI_lower = s$ci_lower,
      CI_upper = s$ci_upper
    )
  })
  
  list(
    boot_distributions = boot_list,   # named list: trait -> numeric vector length B
    summary            = summary_df   # one row per trait
  )
}

#' Bootstrap H² — per time point (matches calc_H2_timeseries)
#'
#' @param data         data.frame
#' @param trait        single trait name (character)
#' @param time_col     name of the time/day column (default "Day")
#' @param time_groups  optional vector of time points; defaults to all unique values
#' @param B            number of bootstrap replicates (default 1000)
#' @param seed         random seed (applied per day for reproducibility; default 123)
#' @return tibble with columns: Time, Trait, H2_mean, H2_se, CI_lower, CI_upper
bootstrap_H2_timeseries <- function(data, trait, time_col = "Day",
                                    time_groups = NULL, B = 100, seed = 123) {
  time_groups <- time_groups %||% sort(unique(data[[time_col]]))
  
  purrr::map_dfr(time_groups, function(d) {
    message("Bootstrapping Day: ", d)
    set.seed(seed)
    df_day  <- dplyr::filter(data, .data[[time_col]] == d)
    h2_boot <- replicate(B, .boot_one(
      data        = df_day,
      trait       = trait,
      formula_str = "~ (1|Genotype)",
      vc_names    = c("Genotype", "Residual")
    ))
    s <- .summarise_boot(h2_boot)
    tibble(
      Time     = d,
      Trait    = trait,
      H2_mean  = s$mean,
      H2_se    = s$se,
      CI_lower = s$ci_lower,
      CI_upper = s$ci_upper
    )
  })
}





################## older version
# #---------------- calc_H2_1dim --------------------------#
# calc_H2 <- function(data, traits) {
#   # the model is use Genotype and Genotype:Replicate as random effect
#   results <- list()
#   
#   for (trait in traits) {
#     # Fit model
#     model <- lmer(
#       as.formula(paste(trait, "~ (1|Genotype)")),
#       data = data
#     )
#     
#     # Extract variance components
#     var_comp <- as.data.frame(VarCorr(model)) %>%
#       dplyr::select(grp, vcov) %>%
#       tidyr::pivot_wider(names_from = grp, values_from = vcov)
#     
#     # Ensure missing components become 0
#     if (!"Genotype" %in% colnames(var_comp)) var_comp$Genotype <- 0
#     if (!"Residual" %in% colnames(var_comp)) var_comp$Residual <- 0
#     
#     # Variance components
#     sigma_G2 <- var_comp$Genotype
#     sigma_E2 <- var_comp$Residual
#     
#     # Broad-sense H²
#     h2 <- round(sigma_G2 / (sigma_G2 + sigma_E2), 4)
#     
#     # Store results
#     results[[trait]] <- data.frame(
#       Trait = trait,
#       Var_Genotype = sigma_G2,
#       Var_Residual = var_comp$Residual,
#       H2 = h2
#     )
#   }
#   
#   # Combine all results into one dataframe
#   h2_results <- do.call(rbind, results)
#   rownames(h2_results) <- NULL
#   
#   return(h2_results)
# }
# 
# 
# 
# # -------------------------- calc_H2_fixedtime --------------------------#
# calc_H2_fixedtime <- function(data, traits, day_col = "Day") {
# # the model is use Genotype and Genotype:Replicate as random effect
#   results <- list()
#   
#   for (trait in traits) {
#     message("Processing: ", trait)
#     
#     # Fit model
#     model <- lmer(
#       as.formula(paste(trait, "~", day_col, "+ (1|Genotype) + (1|Genotype:Replicate)")),
#       data = data
#     )
#     
#     # Extract variance components
#     var_comp <- as.data.frame(VarCorr(model)) %>%
#       dplyr::select(grp, vcov) %>%
#       tidyr::pivot_wider(names_from = grp, values_from = vcov)
#     
#     # Ensure missing components become 0
#     if (!"Genotype" %in% colnames(var_comp)) var_comp$Genotype <- 0
#     if (!"Genotype:Replicate" %in% colnames(var_comp)) var_comp$`Genotype:Replicate` <- 0
#     if (!"Residual" %in% colnames(var_comp)) var_comp$Residual <- 0
#     
#     # Variance components
#     sigma_G2 <- var_comp$Genotype
#     sigma_E2 <- var_comp$`Genotype:Replicate` + var_comp$Residual
#     
#     # Broad-sense H²
#     h2 <- round(sigma_G2 / (sigma_G2 + sigma_E2), 4)
#     
#     # Store results
#     results[[trait]] <- data.frame(
#       Trait = trait,
#       Var_Genotype = sigma_G2,
#       Var_Replicate = var_comp$`Genotype:Replicate`,
#       Var_Residual = var_comp$Residual,
#       H2 = h2
#     )
#   }
#   
#   # Combine all results into one dataframe
#   h2_results <- do.call(rbind, results)
#   rownames(h2_results) <- NULL
#   
#   return(h2_results)
# }
# 
# #-------------------------- calc_H2_timeseries --------------------------#
# calc_H2_timeseries <- function(data, traits, time_col = "Day", time_groups = NULL) {
#   if (is.null(time_groups)) {
#     time_groups <- sort(unique(data[[time_col]]))
#   }
#   
#   results_timeseries <- list()
#   i <- 0
#   
#   for (d in time_groups) {
#     i <- i + 1
#     message("Processing Day: ", d)
#     
#     df_day <- subset(data, data[[time_col]] == d)
#     day_results <- list()
#     
#     for (trait in traits) {
#       cat("  Trait:", trait, "\n")
#       
#       # Fit LMM: (1|Genotype)
#       model <- tryCatch(
#         lmer(as.formula(paste(trait, "~ (1|Genotype)")), data = df_day),
#         error = function(e) {
#           warning("Model failed for trait ", trait, " on day ", d)
#           return(NULL)
#         }
#       )
#       
#       if (is.null(model)) next
#       
#       # Extract variance components
#       var_comp <- as.data.frame(VarCorr(model)) %>%
#         dplyr::select(grp, vcov) %>%
#         tidyr::pivot_wider(names_from = grp, values_from = vcov)
#       
#       # Ensure missing variance terms are 0
#       if (!"Genotype" %in% names(var_comp)) var_comp$Genotype <- 0
#       if (!"Residual" %in% names(var_comp)) var_comp$Residual <- 0
#       
#       # Variance components
#       sigma_G2 <- var_comp$Genotype
#       sigma_E2 <- var_comp$Residual
#       
#       # Broad-sense heritability
#       h2 <- round(sigma_G2 / (sigma_G2 + sigma_E2), 4)
#       
#       # Store result for this trait and day
#       day_results[[trait]] <- data.frame(
#         Time = d,
#         Trait = trait,
#         Var_Genotype = sigma_G2,
#         Var_Residual = sigma_E2,
#         H2 = h2
#       )
#     }
#     
#     # Combine results for this date
#     results_timeseries[[i]] <- do.call(rbind, day_results)
#   }
#   
#   # Combine all days
#   h2_df_timeseries <- do.call(rbind, results_timeseries)
#   rownames(h2_df_timeseries) <- NULL
#   
#   return(h2_df_timeseries)
# }

