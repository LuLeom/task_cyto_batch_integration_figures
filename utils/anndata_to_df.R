#source("helper_functions.R")

#' This function is used to convert an (integrated censored) anndata object
#' into a data.frame() suitable for ggplot plotting.
#' Metadata information is fetched from the full unintegrated dataset
#'
#' @param i_adata AnnData object, integrated data
#' @param u_adata AnnData object, unintegrated dataset
#' @param split_id numeric, split id of the integrated data
#' @return A data.frame() object. The last 7 columns are metadata information, the rest are marker expression.
anndata_to_df <- function(i_adata, u_adata, split_id){
  
  out_adata <- get_obs_var_for_integrated(i_adata, u_adata, split_id)
  out_adata <- subset_markers_tocorrect(out_adata)
  
  out_df <- data.frame(out_adata$layers['integrated'])
  colnames(out_df) <- paste0(out_adata$var[['marker']]," (",out_adata$var[['channel']],")")
  out_df <- cbind(out_df,out_adata$obs)
}
