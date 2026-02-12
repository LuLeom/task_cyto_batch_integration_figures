library(dplyr)
library(tidyr)
library(ComplexHeatmap)

EMDHor_heatmap <- function(df){
# Prep data ####################################################################
marker_cols <- dplyr::setdiff(colnames(df), c("cell_type", "donor"))

mat <- df %>%
  select(all_of(marker_cols)) %>%
  as.matrix()

# Create row labels
row_labels <- paste(df$donor, df$cell_type, sep = "_")

rownames(mat) <- row_labels


# Define row gaps by donor #####################################################
donor_sizes <- df %>%
  count(donor) %>%
  pull(n)

# Gaps after each donor group (except last)
row_gaps <- cumsum(donor_sizes)
row_gaps <- row_gaps[-length(row_gaps)]

# Draw heatmap #################################################################
p <- ComplexHeatmap::Heatmap(mat,
             name = "EMD horizontal",
             col = rev(viridis::magma(100)),
             row_split = df$donor,
             row_gap = unit(1, "mm"),
             show_row_names = TRUE,
             row_names_gp = grid::gpar(fontsize = 6),
             column_names_gp = grid::gpar(fontsize = 6),
             column_title = "Markers",
             row_title = "Donor_cellType",
             heatmap_legend_param = list(title = "EMD horizontal"),
             cluster_rows = FALSE, cluster_columns = FALSE)
return(p)
}

EMDVert_heatmap <- function(df){
  # Prep data ####################################################################
  marker_cols <- dplyr::setdiff(colnames(df), c("cell_type", 
                                                "first_sample", "second_sample"))
  
  mat <- df %>%
    select(all_of(marker_cols)) %>%
    as.matrix()
  
  # Create row labels
  row_labels <- paste(df$cell_type, df$first_sample, df$second_sample, sep = "_")
  
  rownames(mat) <- row_labels
  
  
  # Define row gaps by donor #####################################################
  cellType_sizes <- df %>%
    count(cell_type) %>%
    pull(n)
  
  # Gaps after each donor group (except last)
  row_gaps <- cumsum(cellType_sizes)
  row_gaps <- row_gaps[-length(row_gaps)]
  
  # Draw heatmap #################################################################
  p <- ComplexHeatmap::Heatmap(mat,
                               name = "EMD vertical",
                               col = rev(viridis::magma(100)),
                               row_split = df$cell_type,
                               row_gap = unit(1, "mm"),
                               show_row_names = TRUE,
                               row_names_gp = grid::gpar(fontsize = 6),
                               column_names_gp = grid::gpar(fontsize = 6),
                               column_title = "Markers",
                               row_title = "cellType_first_secondSample",
                               heatmap_legend_param = list(title = "EMD vertical"),
                               cluster_rows = FALSE, cluster_columns = FALSE)
  return(p)
}