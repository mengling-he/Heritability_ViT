
plot_surface_z <- function(
    df_summary,
    z_value  = "neg_var_h2",   # column used for Z-axis (surface height)
    mean_col = "mean_h2",      # column used for surface colors
    highlight_coords = data.frame(
      n_geno = c(20, 30, 40, 60),
      n_rep  = c(6,  4,  3,  2)
    ),
    outdir=NA,
    file_name = "surface_plot.html"
) {
  
  # Axes
  x_vals  <- sort(unique(df_summary$n_geno))
  y_vals  <- sort(unique(df_summary$n_rep))
  metrics <- sort(unique(df_summary$metric))
  
  # Storage
  z_mats    <- list()
  mean_mats <- list()
  highlight_list <- list()
  
  # --- Build matrices for each metric ---
  for (m in metrics) {
    m_df <- df_summary %>% dplyr::filter(metric == m)
    
    z_mat    <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
    mean_mat <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
    
    for (i in seq_along(x_vals)) {
      for (j in seq_along(y_vals)) {
        sub_df <- m_df %>%
          dplyr::filter(n_geno == x_vals[i], n_rep == y_vals[j])
        if (nrow(sub_df) == 1) {
          z_mat[j, i]    <- sub_df[[z_value]]
          mean_mat[j, i] <- sub_df[[mean_col]]
        }
      }
    }
    
    z_mats[[m]]    <- z_mat
    mean_mats[[m]] <- mean_mat
    
    # highlight points for this metric
    highlight_list[[m]] <- highlight_coords %>%
      dplyr::left_join(m_df, by = c("n_geno", "n_rep"))
  }
  
  # --- Plot ---
  p <- plotly::plot_ly()
  n_metrics <- length(metrics)
  
  for (k in seq_along(metrics)) {
    m <- metrics[k]
    z_mat    <- z_mats[[m]]
    mean_mat <- mean_mats[[m]]
    hl       <- highlight_list[[m]]
    
    visible_init <- (k == 1)
    
    # Surface (height = z_value, color = mean_col)
    p <- p %>%
      plotly::add_surface(
        x = x_vals,
        y = y_vals,
        z = z_mat,
        surfacecolor = mean_mat,
        colorscale = "Viridis",
        showscale = (k == 1),
        colorbar = list(title = mean_col),
        visible = visible_init,
        name = paste(m, "surface"),
        hovertemplate = paste(
          "n_geno: %{x}<br>",
          "n_rep: %{y}<br>",
          sprintf("%s: %%{z}", z_value),
          "<extra></extra>"
        )
      )
    
    # Highlight path (line + markers)
    p <- p %>%
      plotly::add_trace(
        data = hl,
        x = ~n_geno,
        y = ~n_rep,
        z = hl[[z_value]],
        type = "scatter3d",
        mode = "lines+markers",
        marker = list(size = 6, color = "red"),
        line   = list(color = "red", width = 5),
        visible = visible_init,
        inherit = FALSE,
        showlegend = FALSE,
        hovertemplate = paste(
          "n_geno: %{x}<br>",
          "n_rep: %{y}<br>",
          sprintf("%s: %%{z}", z_value),
          "<extra></extra>"
        )
      )
    
    # Text labels using mean_col
    p <- p %>%
      plotly::add_text(
        data = hl,
        x = ~n_geno,
        y = ~n_rep,
        z = hl[[z_value]],
        text = sprintf("%s=%.3f", mean_col, hl[[mean_col]]),
        type = "scatter3d",
        mode = "text",
        visible = visible_init,
        inherit = FALSE,
        showlegend = FALSE
      )
  }
  
  # --- Dropdown buttons ---
  n_traces <- 3 * n_metrics
  buttons <- lapply(seq_len(n_metrics), function(k) {
    vis <- rep(FALSE, n_traces)
    idx <- ((k - 1) * 3 + 1):(k * 3)
    vis[idx] <- TRUE
    list(
      method = "update",
      args = list(
        list(visible = vis),
        list(title = paste("Surface of", z_value, "-", metrics[k]))
      ),
      label = metrics[k]
    )
  })
  
  # Layout
  p <- p %>%
    plotly::layout(
      title = paste("Surface of", z_value, "-", metrics[1]),
      scene = list(
        xaxis = list(title = "n_geno"),
        yaxis = list(title = "n_rep"),
        zaxis = list(title = z_value)
      ),
      updatemenus = list(
        list(
          type = "dropdown",
          x = 0.05,
          y = 1.05,
          showactive = TRUE,
          buttons = buttons
        )
      )
    )
  
  # Save only if outdir is provided
  if (!is.na(outdir)) {
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }
    htmlwidgets::saveWidget(p, file.path(outdir, file_name))
  }
  
  return(p)
}







plot_heatmap <- function(df,
                         index_col = "Trait",
                         columns_col = "Time",
                         values_col = "H2",
                         title = "Heritability (H²) by Day and Metric",
                         figsize = c(10, 6),
                         cmap = "YlGnBu",
                         annot = TRUE,
                         fmt = "%.2f",
                         cbar_label = "H²",
                         show = TRUE,
                         save_path = NULL) {
  
  # Pivot data to long format if not already
  df_melt <- df[, c(index_col, columns_col, values_col)]
  colnames(df_melt) <- c("Trait", "Time", "H2")
  
  # Ensure factors are ordered by appearance
  df_melt$Trait <- factor(df_melt$Trait, levels = unique(df_melt$Trait))
  df_melt$Time <- factor(df_melt$Time, levels = unique(df_melt$Time))
  
  # Create heatmap
  p <- ggplot(df_melt, aes(x = Time, y = Trait, fill = H2)) +
    geom_tile(color = "gray80", linewidth = 0.4) +
    scale_fill_distiller(palette = cmap, direction = 1, name = cbar_label) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
    ) +
    labs(title = title, x = columns_col, y = index_col)
  
  # Add annotations if requested
  if (annot) {
    p <- p + geom_text(aes(label = sprintf(fmt, H2)), size = 3, color = "black")
  }
  
  # Adjust figure size and show/save
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = figsize[1], height = figsize[2])
  }
  
  if (show) print(p)
  
  return(p)
}