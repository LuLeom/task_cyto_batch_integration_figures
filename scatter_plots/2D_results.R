library(ggplot2)
library(ggrepel)

data <- read.csv("datasets/260223_cytobenchmark_run/funky_heatmap/overall_scores.csv")

control_methods <- c( "shuffle_integration_globally",
                     "shuffle_integration_within_cell_type", "shuffle_integration_within_batch")      

data <- data[!data$method %in% control_methods,]
ggplot(data, aes(x = mean_batch_mixing, y = mean_bio_conservation, label = method)) +
  geom_point() +
  geom_text_repel(max.overlaps = 100) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Batch mixing", y = "Bio-conservation") +
  theme_bw() +
  theme(panel.grid = element_blank())
