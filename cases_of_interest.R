library(ggplot2)
library(dplyr)
library(anndata)
source("utils/anndata_to_df.R")
source("utils/helper_functions.R")
source("density_plots/EMD_densities.R")
source("heatmap_plots/EMD_heatmaps.R")


methods <- c("cytonorm_no_controls_to_goal", "cytonorm_no_controls_to_mid")
data_path <- "datasets/260223_cytobenchmark_run/intermediate_files/mouse_spleen_flow_cytometry/method_out/" # Change with your local path
                  
##### Neutrophils CD64 MSFC
EMDHor_densities("cytonorm_no_controls_to_goal", 
                 cols_to_plot = c("CD64#BV711 (V710-A)"), 
                 cellType = "Neutrophils",
                 integrated_dataset_s1 = paste0(data_path,"cytonorm_no_controls_to_goal_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cytonorm_no_controls_to_goal_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/Neutrophils_CD64_CN_noCtrls_toGoal.pdf")

EMDHor_densities("cycombine_no_controls_to_goal", 
                 cols_to_plot = c("CD64#BV711 (V710-A)"), 
                 cellType = "Neutrophils",
                 integrated_dataset_s1 = paste0(data_path,"cycombine_no_controls_to_goal_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cycombine_no_controls_to_goal_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/Neutrophils_CD64_CC_noCtrls_toGoal.pdf")

EMDHor_densities("cytonorm_no_controls_to_mid", 
                 cols_to_plot = c("CD64#BV711 (V710-A)"), 
                 cellType = "Neutrophils",
                 integrated_dataset_s1 = paste0(data_path,"cytonorm_no_controls_to_mid_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cytonorm_no_controls_to_mid_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/Neutrophils_CD64_CN_noCtrls_toMid.pdf")

EMDHor_densities("cycombine_no_controls_to_mid", 
                 cols_to_plot = c("CD64#BV711 (V710-A)"), 
                 cellType = "Neutrophils",
                 integrated_dataset_s1 = paste0(data_path,"cycombine_no_controls_to_mid_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cycombine_no_controls_to_mid_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/Neutrophils_CD64_CC_noCtrls_toMid.pdf")

##### Neutrophils CD16 MSFC 
EMDHor_densities("cytonorm_no_controls_to_goal", 
                 cols_to_plot = c("CD161#BV605 (V605-A)"), 
                 cellType = "Neutrophils",
                 integrated_dataset_s1 = paste0(data_path,"cytonorm_no_controls_to_goal_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cytonorm_no_controls_to_goal_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/Neutrophils_CD161_CN_noCtrls_toGoal.pdf")

EMDHor_densities("cycombine_no_controls_to_goal", 
                 cols_to_plot = c("CD161#BV605 (V605-A)"), 
                 cellType = "Neutrophils",
                 integrated_dataset_s1 = paste0(data_path,"cycombine_no_controls_to_goal_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cycombine_no_controls_to_goal_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/Neutrophils_CD161_CC_noCtrls_toGoal.pdf")

EMDHor_densities("cytonorm_no_controls_to_mid", 
                 cols_to_plot = c("CD161#BV605 (V605-A)"), 
                 cellType = "Neutrophils",
                 integrated_dataset_s1 = paste0(data_path,"cytonorm_no_controls_to_mid_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cytonorm_no_controls_to_mid_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/Neutrophils_CD161_CN_noCtrls_toMid.pdf")

EMDHor_densities("cycombine_no_controls_to_mid", 
                 cols_to_plot = c("CD161#BV605 (V605-A)"), 
                 cellType = "Neutrophils",
                 integrated_dataset_s1 = paste0(data_path,"cycombine_no_controls_to_mid_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cycombine_no_controls_to_mid_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/Neutrophils_CD161_CC_noCtrls_toMid.pdf")

##### CD4 marker fro CD4 Naive cells / "CD4#BUV496 (UV515-A)"
EMDHor_densities("cytonorm_no_controls_to_goal", 
                 cols_to_plot = c("CD4#BUV496 (UV515-A)"), 
                 cellType = "CD4 Naive",
                 integrated_dataset_s1 = paste0(data_path,"cytonorm_no_controls_to_goal_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cytonorm_no_controls_to_goal_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/CD4Naive_CD4_CN_noCtrls_toGoal.pdf")

EMDHor_densities("cytonorm_no_controls_to_mid", 
                 cols_to_plot = c("CD4#BUV496 (UV515-A)"), 
                 cellType = "CD4 Naive",
                 integrated_dataset_s1 = paste0(data_path,"cytonorm_no_controls_to_mid_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cytonorm_no_controls_to_mid_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/CD4Naive_CD4_CN_noCtrls_toMid.pdf")

EMDHor_densities("cycombine_no_controls_to_goal", 
                 cols_to_plot = c("CD4#BUV496 (UV515-A)"), 
                 cellType = "CD4 Naive",
                 integrated_dataset_s1 = paste0(data_path,"cycombine_no_controls_to_goal_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cycombine_no_controls_to_goal_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/CD4Naive_CD4_CC_noCtrls_toGoal.pdf")

EMDHor_densities("cycombine_no_controls_to_mid", 
                 cols_to_plot = c("CD4#BUV496 (UV515-A)"), 
                 cellType = "CD4 Naive",
                 integrated_dataset_s1 = paste0(data_path,"cycombine_no_controls_to_mid_split1.h5ad"),
                 integrated_dataset_s2 = paste0(data_path,"cycombine_no_controls_to_mid_split2.h5ad"),
                 unintegrated_dataset = "datasets/unintegrated/mouse_spleen_flow_cytometry.h5ad",
                 fileName = "output/CD4Naive_CD4_CC_noCtrls_toMid.pdf")
