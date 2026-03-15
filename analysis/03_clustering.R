# 03_clustering.R
# Find neighbors, cluster at multiple resolutions, generate UMAPs

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)

out_dir <- "data/processed"
dir.create("results/clustering", showWarnings = FALSE, recursive = TRUE)

# --- Load ---
cat("Loading merged object...\n")
merged <- readRDS(file.path(out_dir, "merged.rds"))

# --- Elbow plot ---
p_elbow <- ElbowPlot(merged, ndims = 50) + ggtitle("Elbow Plot")
ggsave("results/clustering/elbow_plot.pdf", p_elbow, width = 7, height = 5)

# --- Clustering ---
cat("Finding neighbors...\n")
merged <- FindNeighbors(merged, dims = 1:30, verbose = FALSE)

cat("Clustering at multiple resolutions...\n")
merged <- FindClusters(merged, resolution = 0.3, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.8, verbose = FALSE)

# Default to resolution 0.5
Idents(merged) <- "SCT_snn_res.0.5"

# --- UMAP plots at each resolution ---
p_03 <- DimPlot(merged, group.by = "SCT_snn_res.0.3", label = TRUE, repel = TRUE) +
  ggtitle("Resolution 0.3") + NoLegend()
p_05 <- DimPlot(merged, group.by = "SCT_snn_res.0.5", label = TRUE, repel = TRUE) +
  ggtitle("Resolution 0.5") + NoLegend()
p_08 <- DimPlot(merged, group.by = "SCT_snn_res.0.8", label = TRUE, repel = TRUE) +
  ggtitle("Resolution 0.8") + NoLegend()

ggsave("results/clustering/umap_resolutions.pdf", p_03 + p_05 + p_08,
       width = 18, height = 6)

# Cluster counts
cat("\nCluster sizes (res 0.5):\n")
print(table(merged$SCT_snn_res.0.5))

cat("\nCluster sizes by condition (res 0.5):\n")
print(table(merged$SCT_snn_res.0.5, merged$condition))

# --- Save ---
saveRDS(merged, file.path(out_dir, "merged_clustered.rds"))

cat("Clustering complete.\n")
