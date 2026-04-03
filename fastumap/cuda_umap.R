library(ggplot2)

plot_umap <- function(file, title) {
  dat <- read.csv(file)
  ggplot(dat, aes(x = umap_x, y = umap_y)) +
    geom_point(size = 0.3, alpha = 0.5) +
    theme_minimal() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2")
}

plot_umap("fastumap/umap_mock5.csv", "mock5 - cosine")
plot_umap("fastumap/umap_mock5_sign.csv", "mock5 - cosine + sign")
plot_umap("fastumap/umap_mock5_log.csv", "mock5 - cosine + log")
plot_umap("fastumap/umap_mock5_hamming.csv", "mock5 - hamming + sign")
plot_umap("fastumap/umap_mock5_n30.csv", "mock5 - cosine + n_neighbors=30")

plot_umap("fastumap/umap_s1.csv", "s1 - cosine")
plot_umap("fastumap/umap_s1_sign.csv", "s1 - cosine + sign")
plot_umap("fastumap/umap_s1_log.csv", "s1 - cosine + log")
plot_umap("fastumap/umap_s1_hamming.csv", "s1 - hamming + sign")
plot_umap("fastumap/umap_s1_n30.csv", "s1 - cosine + n_neighbors=30")

plot_umap("fastumap/umap_s2.csv", "s2 - cosine")
plot_umap("fastumap/umap_s2_sign.csv", "s2 - cosine + sign")
plot_umap("fastumap/umap_s2_log.csv", "s2 - cosine + log")
plot_umap("fastumap/umap_s2_hamming.csv", "s2 - hamming + sign")
plot_umap("fastumap/umap_s2_n30.csv", "s2 - cosine + n_neighbors=30")

plot_umap("fastumap/umap_s2_n100.csv", "s2 - cosine + n_neighbors=100")
plot_umap("fastumap/umap_s2_n500.csv", "s2 - cosine + n_neighbors=500")
plot_umap("fastumap/umap_s2_n1000.csv", "s2 - cosine + n_neighbors=1000")

plot_umap("fastumap/umap_s2_sign_n1000.csv", "s2 - cosine + sign + n_neighbors=1000")
plot_umap("fastumap/umap_s2_sign2_n1000.csv", "s2 - cosine + sign2 + n_neighbors=1000") + xlim(35,50) + ylim(110,120)
