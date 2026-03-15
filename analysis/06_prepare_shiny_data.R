# 06_prepare_shiny_data.R
# Export lightweight data for the Shiny app (no Seurat dependency needed)

library(Seurat)
library(Matrix)
library(dplyr)
library(RColorBrewer)

set.seed(42)

proc_dir <- "data/processed"
app_dir <- "app/data"
dir.create(app_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load ---
cat("Loading annotated object...\n")
merged <- readRDS(file.path(proc_dir, "merged_annotated.rds"))

# --- 1. UMAP + metadata ---
cat("Exporting UMAP + metadata...\n")
umap_coords <- Embeddings(merged, "umap")
meta <- merged@meta.data

umap_metadata <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  cell_type = meta$cell_type,
  condition = meta$condition,
  cluster = as.character(Idents(merged)),
  nFeature_RNA = meta$nFeature_RNA,
  nCount_RNA = meta$nCount_RNA,
  percent_mt = meta$percent.mt,
  stringsAsFactors = FALSE
)
saveRDS(umap_metadata, file.path(app_dir, "umap_metadata.rds"))

# --- 2. Expression matrix (SCT normalized, sparse) ---
cat("Exporting expression matrix...\n")
expr_matrix <- GetAssayData(merged, assay = "SCT", layer = "data")
# Keep as sparse matrix
saveRDS(expr_matrix, file.path(app_dir, "expression_matrix.rds"))
cat(sprintf("Expression matrix: %d genes x %d cells (%.1f MB sparse)\n",
            nrow(expr_matrix), ncol(expr_matrix),
            object.size(expr_matrix) / 1e6))

# --- 3. DE results ---
cat("Exporting DE results...\n")
de_results <- readRDS(file.path(proc_dir, "de_results.rds"))
saveRDS(de_results, file.path(app_dir, "de_results.rds"))

# --- 4. Cluster markers ---
cat("Exporting marker genes...\n")
marker_genes <- readRDS(file.path(proc_dir, "all_cluster_markers.rds"))
saveRDS(marker_genes, file.path(app_dir, "marker_genes.rds"))

# --- 5. Gene list ---
cat("Exporting gene list...\n")
gene_list <- sort(rownames(expr_matrix))
saveRDS(gene_list, file.path(app_dir, "gene_list.rds"))
cat(sprintf("Gene list: %d genes\n", length(gene_list)))

# --- 6. Color palette ---
cat("Creating color palette...\n")
cell_types <- sort(unique(umap_metadata$cell_type))
n_types <- length(cell_types)

# Use a combination of palettes for distinct colors
if (n_types <= 12) {
  colors <- brewer.pal(max(n_types, 3), "Set3")[1:n_types]
} else {
  colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_types)
}

# Manual color assignment for better aesthetics
color_map <- c(
  "Type IIB myonuclei"    = "#E41A1C",
  "Type IIX myonuclei"    = "#FF7F00",
  "Type IIA myonuclei"    = "#FFC107",
  "Type I myonuclei"      = "#984EA3",
  "Denervated myonuclei"  = "#1B9E77",
  "Catabolic myonuclei"   = "#D95F02",
  "Satellite cells"       = "#4DAF4A",
  "Endothelial"           = "#377EB8",
  "Lymphatic endothelial" = "#6BAED6",
  "Fibroblasts/FAPs"      = "#A65628",
  "Tenocytes"             = "#C49A6C",
  "Immune"                = "#F781BF",
  "Pericytes"             = "#999999",
  "Adipocytes"            = "#66C2A5",
  "NMJ/MTJ nuclei"       = "#8DD3C7",
  "Unknown"               = "#CCCCCC"
)

# Use manual colors where available, fill in rest
color_palette <- setNames(colors, cell_types)
for (ct in names(color_map)) {
  if (ct %in% cell_types) {
    color_palette[ct] <- color_map[ct]
  }
}

saveRDS(color_palette, file.path(app_dir, "color_palette.rds"))

# --- 7. Summary stats ---
cat("Exporting summary stats...\n")
summary_stats <- umap_metadata %>%
  group_by(cell_type, condition) %>%
  summarise(
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_umis = median(nCount_RNA),
    .groups = "drop"
  )

saveRDS(summary_stats, file.path(app_dir, "summary_stats.rds"))

# --- Report ---
cat("\n=== Shiny data export summary ===\n")
cat(sprintf("Total cells: %d\n", nrow(umap_metadata)))
cat(sprintf("Total genes: %d\n", length(gene_list)))
cat(sprintf("Cell types: %s\n", paste(cell_types, collapse = ", ")))
cat(sprintf("DE cell types: %s\n", paste(names(de_results), collapse = ", ")))

files <- list.files(app_dir, full.names = TRUE)
for (f in files) {
  cat(sprintf("  %s: %.1f MB\n", basename(f), file.info(f)$size / 1e6))
}

cat("\nShiny data preparation complete.\n")
