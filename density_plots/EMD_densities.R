library(ggplot2)
library(dplyr)
library(anndata)

EMDVert_densities <- function(methods,
                              cols_to_plot, 
                              cellType = "all",
                              integrated_dataset_s1,
                              integrated_dataset_s2,
                              unintegrated_dataset,
                              fileName,
                              manXlim = c(-1,5), 
                              show_controls = FALSE,
                              show_group = TRUE){
  pList <- EMDVert_makePlots(methods = methods, 
                             cols_to_plot = cols_to_plot, 
                             cellType = cellType,
                             integrated_dataset_s1 = integrated_dataset_s1,
                             integrated_dataset_s2 = integrated_dataset_s2,
                             unintegrated_dataset = unintegrated_dataset,
                             manXlim = manXlim, 
                             show_controls = show_controls,
                             show_group = show_group)
  
  method_splits <- names(pList)
  markers <- names(pList[[1]])
  
  plot_grid_list <- unlist(
    lapply(method_splits, function(ms) {
      lapply(markers, function(marker) {
        pList[[ms]][[marker]]
      })
    }),
    recursive = FALSE
  )
  
  pdf(fileName, height = length(method_splits)*3, width = length(markers)*3)
  print(ggpubr::ggarrange(
    plotlist = plot_grid_list,
    nrow = length(method_splits),
    ncol = length(markers),
    common.legend = TRUE,
    legend = "left"))
  dev.off()  
}

EMDVert_makePlots <- function(methods, 
                              cols_to_plot, 
                              cellType, 
                              integrated_dataset_s1,
                              integrated_dataset_s2,
                              unintegrated_dataset,
                              manXlim, 
                              show_controls, 
                              show_group){
  pList <- list()
  
  # Make plot for unintegrated data
  unintegrated <- anndata::read_h5ad(unintegrated_dataset)
  
  df <- unintegrated$layers['preprocessed']
  colnames(df) <- paste0(unintegrated$var[["marker"]], " (", unintegrated$var[["channel"]], ")")
  tmp <- densityPlots(df = data.frame(cbind(df,
                                            unintegrated$obs),
                                      check.names = FALSE),
                      cols_to_plot = cols_to_plot, 
                      cellType = cellType,
                      show_controls = show_controls,
                      show_group = show_group,
                      manXlim = manXlim)
  for (i in 1:length(tmp)){
    tmp[[i]] <- tmp[[i]] + ggtitle(cols_to_plot[i])
  }
  tmp[[1]] <- tmp[[1]] + ylab("unintegrated")
  
  pList[["unintegrated"]] <- tmp
  
  # Make plots for methods
  for (method in methods){
    writeLines(method)
    split <- 0 #temp 0 as placeholder
    for (dataset in c(integrated_dataset_s1, integrated_dataset_s2)){
      split <- split + 1
      writeLines(as.character(split))
      df <- anndata_to_df(i_adata = anndata::read_h5ad(dataset),
                          u_adata = unintegrated,
                          split_id = split)
      df$group <- factor(df$group)
      method_short <- sub("batchadjust", "BA", 
                          sub("cycombine", "CC", 
                              sub("cytonorm", "CN",
                                  sub("_control", "Ctrl", 
                                      sub("_controls", "Ctrls", method)))))
      method_short <- paste0(method_short, "_", sub("split", "", split))
      tmp <- densityPlots(df, 
                          cols_to_plot = cols_to_plot, 
                          cellType = cellType,
                          show_controls = show_controls,
                          show_group = show_group,
                          manXlim = manXlim)
      tmp[[1]] <- tmp[[1]] + ylab(method_short)
      pList[[method_short]] <- tmp
    }
  } 
  return(pList)
}

densityPlots <- function(df, 
                         cols_to_plot, 
                         cellType,
                         show_controls, 
                         show_group, 
                         manXlim) {
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
    if (show_group){
      ggplot(df, aes(x = .data[[marker_col]])) +
        geom_density(aes(color = batch, 
                         group = interaction(sample, batch),
                         linetype = group)) +
        theme_minimal() +
        xlim(manXlim) +
        xlab("") + ylab("") +
        scale_linetype(drop = FALSE)
    } else {
      ggplot(df, aes(x = .data[[marker_col]])) +
        geom_density(aes(color = batch, 
                         group = interaction(sample, batch))) +
        theme_minimal() +
        xlim(manXlim) +
        xlab("") + ylab("")
    }
  })
  
  # Name the list by marker
  names(plots) <- marker_cols
  
  return(plots)
}



EMDHor_densities <- function(method,
                             cols_to_plot, 
                             cellType = "all",
                             integrated_dataset_s1,
                             integrated_dataset_s2,
                             unintegrated_dataset,
                             fileName,
                             manXlim = c(-1,5), 
                             show_controls = FALSE,
                             show_group = TRUE){
  pList <- EMDHor_makePlots(method = method, 
                            cols_to_plot = cols_to_plot, 
                            cellType = cellType,
                            integrated_dataset_s1 = integrated_dataset_s1,
                            integrated_dataset_s2 = integrated_dataset_s2,
                            unintegrated_dataset = unintegrated_dataset,
                            manXlim = manXlim, 
                            show_controls = show_controls,
                            show_group = show_group)
  
  donors <- names(pList)
  markers <- names(pList[[1]])
  
  plot_grid_list <- unlist(
    lapply(donors, function(ms) {
      lapply(markers, function(marker) {
        pList[[ms]][[marker]]
      })
    }),
    recursive = FALSE
  )
  
  pdf(fileName, height = length(donors)*3, width = length(markers)*3)
  print(ggpubr::ggarrange(
    plotlist = plot_grid_list,
    nrow = length(donors),
    ncol = length(markers),
    common.legend = TRUE,
    legend = "left"))
  dev.off()  
}

EMDHor_makePlots <- function(method, 
                             cols_to_plot, 
                             cellType, 
                             integrated_dataset_s1,
                             integrated_dataset_s2,
                             unintegrated_dataset,
                             manXlim, 
                             show_controls, 
                             show_group){
  pList <- list()
  
  # Make plot for unintegrated data
  unintegrated <- anndata::read_h5ad(unintegrated_dataset)
  
  df <- unintegrated$layers['preprocessed']
  colnames(df) <- paste0(unintegrated$var[["marker"]], " (", unintegrated$var[["channel"]], ")")
  tmp <- densityPlots(df = data.frame(cbind(df,
                                            unintegrated$obs),
                                      check.names = FALSE),
                      cols_to_plot = cols_to_plot, 
                      cellType = cellType,
                      show_controls = show_controls,
                      show_group = show_group,
                      manXlim = manXlim)
  for (i in 1:length(tmp)){
    tmp[[i]] <- tmp[[i]] + ggtitle(cols_to_plot[i])
  }
  tmp[[1]] <- tmp[[1]] + ylab("unintegrated")
  
  pList[["unintegrated"]] <- tmp
  
  # Make plots for methods
  df <- rbind(anndata_to_df(i_adata = anndata::read_h5ad(integrated_dataset_s1),
                            u_adata = unintegrated, split_id = "split1"),
              anndata_to_df(i_adata = anndata::read_h5ad(integrated_dataset_s2),
                            u_adata = unintegrated, split_id = "split2"))
  df$group <- factor(df$group)
  print(df)
  for (donor in unique(df$donor)){
    writeLines(as.character(donor))
    method_short <- sub("batchadjust", "BA", 
                        sub("cycombine", "CC", 
                            sub("cytonorm", "CN",
                                sub("_control", "Ctrl", 
                                    sub("_controls", "Ctrls", method)))))
    tmp <- densityPlots(df[df$donor == donor,], 
                        cols_to_plot = cols_to_plot, 
                        cellType = cellType,
                        show_controls = show_controls,
                        show_group = show_group,
                        manXlim = manXlim)
    tmp[[1]] <- tmp[[1]] + ylab(donor)
    pList[[donor]] <- tmp
  } 
  return(pList)
}
