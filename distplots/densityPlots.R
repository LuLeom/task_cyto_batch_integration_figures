library(ggplot2)
library(dplyr)

densityPlots <- function(df, cols_to_plot, cellType = "all", manXlim = c(-1,5),
                         show_controls = FALSE) {
  df$batch <- factor(df$batch)
  
  # Optionally subset by cellType
  if (cellType != "all") {
    df <- df %>% filter(cell_type == !!cellType)
  }
  
  # Optionally remove controls
  if (!show_controls) {
    df <- df %>% filter(is_control == 0)
  }
  
  # Find marker columns
  marker_cols <- cols_to_plot
  
  # Create one plot per marker
  plots <- lapply(marker_cols, function(marker_col) {
    
    ggplot(df, aes(x = .data[[marker_col]])) +
      geom_density(aes(color = batch, 
                       group = interaction(sample, batch))) +
      theme_minimal() +
      xlim(manXlim) +
      xlab("") + ylab("")
  })
  
  # Name the list by marker
  names(plots) <- marker_cols
  
  return(plots)
}