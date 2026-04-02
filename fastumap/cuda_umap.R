library(ggplot2)

dat <- read.csv("fastumap/umap_s2.csv")

ggplot(dat, aes(x = umap_x, y = umap_y)) +
  geom_point(size = 0.3, alpha = 0.5) +
  theme_minimal() +
  labs(title = "UMAP (cosine, CUDA)", x = "UMAP 1", y = "UMAP 2")

#ggsave("/tmp/umap_plot.png", width = 8, height = 7, dpi = 150)





dat <- read.csv("fastumap/umap_mock5.csv")

ggplot(dat, aes(x = umap_x, y = umap_y)) +
  geom_point(size = 0.3, alpha = 0.5) +
  theme_minimal() +
  labs(title = "UMAP (cosine, CUDA)", x = "UMAP 1", y = "UMAP 2")

#ggsave("/tmp/umap_plot.png", width = 8, height = 7, dpi = 150)