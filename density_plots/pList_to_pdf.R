pList_to_pdf <- function(methods,
                         cols_to_plot, 
                         cellType = "all",
                         dataset,
                         fileName,
                         manXlim = c(-1,5), 
                         show_controls = FALSE){
  pList <- makePlots(methods = methods, 
                     cols_to_plot = cols_to_plot, 
                     cellType = cellType,
                     dataset = dataset,
                     manXlim = manXlim, 
                     show_controls = show_controls)
  
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
