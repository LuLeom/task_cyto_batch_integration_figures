source("helper_functions.R")
anndata_to_df <- function(i_adata, u_adata, split_id){
  meta <- c("cell_type","batch","sample","donor","group","is_control","split")
  
  out_adata <- get_obs_var_for_integrated(i_adata, u_adata, split_id)
}