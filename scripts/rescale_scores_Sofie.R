### Functions ##################################################################
select_dataset <- function(df, dataset){
  res <- df[df$dataset == dataset, -c(1,2)]
  rownames(res) <- df[df$dataset == dataset, "method"]
  res
}

scale_scores <- function(scores){
  scores_max <- scores["perfect_integration", ]
  scores_min <- scores["no_integration",]
  scores_min["ratio_consistent_peaks"] <- 0.5
  
  to_invert <- which(scores_max < scores_min)
  scores[,to_invert] <- -scores[,to_invert]
  
  # recompute after inversion
  scores_max <- scores["perfect_integration", ]
  scores_min <- scores["no_integration",]
  scores_min["ratio_consistent_peaks"] <- 0.5
  
  scores_scaled <- do.call(rbind, apply(scores, 1, function(x) (x - scores_min)/(scores_max - scores_min)))
  scores_scaled[scores_scaled > 1] <- 1
  scores_scaled[scores_scaled < 0] <- 0
  scores_scaled
}

### Scale data #################################################################
scores <- read.csv("data/raw_dataset_scores.csv")

human_scores <- select_dataset(scores, "human_blood_mass_cytometry")
lille_scores <- select_dataset(scores, "lille_spectral_flow_cytometry")
mouse_scores <- select_dataset(scores, "mouse_spleen_flow_cytometry")

human_scaled <- scale_scores(human_scores)
lille_scaled <- scale_scores(lille_scores)
mouse_scaled <- scale_scores(mouse_scores)

### Heatmaps####################################################################
pheatmap::pheatmap(mouse_scaled, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   display_numbers = round(mouse_scaled, 2),
                   scale = "none")

plot(mouse_scores$flowsom_mean_mapping_similarity)

tmp <- array(NA, dim = c(nrow(human_scaled), ncol(human_scaled), 3))
tmp[,,1] <- as.matrix(human_scaled)
tmp[,,2] <- as.matrix(mouse_scaled)
tmp[,,3] <- as.matrix(lille_scaled)

tmp_mean <- apply(tmp, c(1,2), mean, na.rm = TRUE)
dimnames(tmp_mean) <- dimnames(human_scaled)


pheatmap::pheatmap(tmp_mean, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   display_numbers = round(tmp_mean, 2),
                   scale = "none")


pheatmap::pheatmap(mouse_scores, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   display_numbers = round(mouse_scores, 2),
                   scale = "none")

mean_batchmixing <- apply(tmp_mean[,c("average_batch_r2_ct", "emd_mean_ct_horiz", "emd_mean_ct_vert", "iLisi")], 1, mean, na.rm =TRUE)
mean_bioconservation <- apply(tmp_mean[,c("cLisi", "flowsom_mean_mapping_similarity", "ratio_consistent_peaks")], 1, mean, na.rm =TRUE)

df_overview <- data.frame(method = names(mean_batchmixing),
                          mean_batchmixing = mean_batchmixing,
                          mean_bioconservation = mean_bioconservation)

### Overview figure ############################################################
library(ggplot2)
library(ggrepel)
ggplot(df_overview, aes(x = mean_batchmixing, y = mean_bioconservation, label = method)) +
  geom_point() +
  geom_text_repel(max.overlaps = 100, col = "lightgrey") + 
  theme_minimal()

### Line plots showing the original and scaled scores ##########################
library(tidyr)
library(dplyr)
human_scores$method  <- rownames(human_scores)
human_scores$flowsom_mean_mapping_similarity <- human_scores$flowsom_mean_mapping_similarity/100
human_scaled$method  <- rownames(human_scaled)
scores_long <- human_scores %>%
  pivot_longer(-method, names_to = "metric", values_to = "value") %>%
  mutate(type = "Before")
scaled_long <- human_scaled %>%
  pivot_longer(-method, names_to = "metric", values_to = "value") %>%
  mutate(type = "After")
plot_data <- bind_rows(scores_long, scaled_long)
plot_data$type <- factor(plot_data$type,
                         levels = c("Before", "After"))
colors <- c("#F8766D", "#ED8141", "#E08B00", "#CF9400", "#BB9D00", "#A3A500", 
            "#85AD00", "#5BB300", "#00B81F", "#00BC59", "#00BF7D", "#00C19C", 
            "#00C0B8", "#00BDD0", "#00B8E5", "#00B0F6", "#00A5FF", "#7997FF", 
            "#AC88FF", "#CF78FF", "#E76BF3", "#F763E0", "#FF61C9", "#FF65AE", 
            "#FF6C90")
names(colors) <- unique(plot_data$method)
colors["no_integration"] <- "#780000"
colors["perfect_integration"] <- "#003049"

ggplot(plot_data, aes(x = type, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, aes(group = type)) +
  geom_point(aes(color = method, group = method), size = 2) +
  geom_line(alpha = 0.4, aes(group = method, color = method)) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_manual(values = colors) +
  theme_bw() 

plot_data_sub <- plot_data[!grepl("shuffle", plot_data$method),]
ggplot(plot_data_sub, aes(x = type, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, aes(group = type)) +
  geom_point(aes(color = method, group = method), size = 2) +
  geom_line(alpha = 0.4, aes(group = method, color = method)) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_manual(values = colors) +
  theme_bw() 
