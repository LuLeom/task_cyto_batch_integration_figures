library(ggplot2)
library(dplyr)
library(anndata)

makePlots <- function(methods, cols_to_plot, cellType = "all", dataset, 
                      manXlim = c(-1,5), show_controls = FALSE){
  pList <- list()
  
  # Make plot for unintegrated data
  unintegrated <- anndata::read_h5ad(paste0("intermediate_files/", dataset, "/unintegrated.h5ad"))
  
  df <- unintegrated$layers['preprocessed']
  colnames(df) <- paste0(unintegrated$var[["marker"]], " (", unintegrated$var[["channel"]], ")")
  tmp <- densityPlots(df = data.frame(cbind(df,
                                            unintegrated$obs),
                                      check.names = FALSE),
                      cols_to_plot = cols_to_plot, 
                      cellType = cellType,
                      show_controls = show_controls,
                      manXlim = manXlim)
  for (i in 1:length(tmp)){
    tmp[[i]] <- tmp[[i]] + ggtitle(cols_to_plot[i])
  }
  tmp[[1]] <- tmp[[1]] + ylab("unintegrated")
  
  pList[["unintegrated"]] <- tmp
  
  # Make plots for methods
  for (method in c(methods)){
    writeLines(method)
    for (split in c("split1", "split2")){
      writeLines(split)
      df <- anndata_to_df(i_adata = anndata::read_h5ad(paste0("intermediate_files/", dataset, "/method_out/", method, "_", split, ".h5ad")),
                          u_adata = unintegrated,
                          split_id = "split1")
      method_short <- sub("batchadjust", "BA", 
                          sub("cycombine", "CC", 
                              sub("cytonorm", "CN",
                                  sub("_control", "Ctrl", 
                                      sub("_controls", "Ctrls", method)))))
      method_short <- paste0(method_short, "_", sub("split", "", split))
      tmp <- densityPlots(df, 
                          cols_to_plot = cols_to_plot, 
                          cellType = cellType,
                          show_controls = show_controls)
      tmp[[1]] <- tmp[[1]] + ylab(method_short)
      pList[[method_short]] <- tmp
    }
  } 
  return(pList)
}