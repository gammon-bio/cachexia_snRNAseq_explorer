# 02_integration.R
# Merge control + KIC, re-normalize jointly, check for batch effects

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)

out_dir <- "data/processed"
dir.create("results/integration", showWarnings = FALSE, recursive = TRUE)

# --- Load ---
cat("Loading SCTransformed objects...\n")
wt <- readRDS(file.path(out_dir, "wt_sct.rds"))
kic <- readRDS(file.path(out_dir, "kic_sct.rds"))

# --- Merge ---
cat("Merging objects...\n")
merged <- merge(wt, y = kic, add.cell.ids = c("Control", "KIC"),
                project = "ZhangCachexia")

cat(sprintf("Merged: %d nuclei\n", ncol(merged)))

# --- Joint SCTransform on merged object ---
cat("Running joint SCTransform...\n")
merged <- SCTransform(merged, vst.flavor = "v2", verbose = FALSE)

# --- PCA + UMAP for initial visualization ---
merged <- RunPCA(merged, npcs = 50, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)

# Check for batch effects
p1 <- DimPlot(merged, group.by = "condition", reduction = "umap") +
  ggtitle("UMAP by Condition (pre-integration)")
p2 <- DimPlot(merged, group.by = "orig.ident", reduction = "umap") +
  ggtitle("UMAP by Sample")

ggsave("results/integration/umap_pre_integration.pdf", p1 + p2, width = 14, height = 6)

# --- Check if Harmony integration is needed ---
# Simple merge often sufficient for 2 samples from same experiment
# If batch effects are visible, uncomment below:
#
# library(harmony)
# cat("Running Harmony integration...\n")
# merged <- RunHarmony(merged, group.by.vars = "condition",
#                      reduction = "pca", assay.use = "SCT",
#                      reduction.save = "harmony")
# merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30,
#                   reduction.name = "umap.harmony")
# p3 <- DimPlot(merged, group.by = "condition", reduction = "umap.harmony") +
#   ggtitle("UMAP by Condition (post-Harmony)")
# ggsave("results/integration/umap_post_harmony.pdf", p3, width = 7, height = 6)

# --- Save ---
saveRDS(merged, file.path(out_dir, "merged.rds"))

cat("Integration complete.\n")
