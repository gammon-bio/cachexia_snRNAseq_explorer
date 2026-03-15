# test_app_data.R
# Verify pre-computed Shiny data files are correct

library(Matrix)

app_data_dir <- "app/data"

cat("=== Shiny App Data Validation ===\n\n")

# --- Check all files exist ---
required_files <- c(
  "umap_metadata.rds",
  "expression_matrix.rds",
  "de_results.rds",
  "marker_genes.rds",
  "gene_list.rds",
  "color_palette.rds",
  "summary_stats.rds"
)

cat("1. Checking required files...\n")
for (f in required_files) {
  path <- file.path(app_data_dir, f)
  if (file.exists(path)) {
    size_mb <- file.info(path)$size / 1e6
    cat(sprintf("  [OK] %s (%.1f MB)\n", f, size_mb))
  } else {
    cat(sprintf("  [MISSING] %s\n", f))
  }
}

# --- Load and validate ---
cat("\n2. Loading data...\n")
umap_meta <- readRDS(file.path(app_data_dir, "umap_metadata.rds"))
expr_matrix <- readRDS(file.path(app_data_dir, "expression_matrix.rds"))
de_results <- readRDS(file.path(app_data_dir, "de_results.rds"))
marker_genes <- readRDS(file.path(app_data_dir, "marker_genes.rds"))
gene_list <- readRDS(file.path(app_data_dir, "gene_list.rds"))
color_palette <- readRDS(file.path(app_data_dir, "color_palette.rds"))
summary_stats <- readRDS(file.path(app_data_dir, "summary_stats.rds"))

# --- UMAP metadata ---
cat("\n3. UMAP metadata:\n")
cat(sprintf("  Rows: %d\n", nrow(umap_meta)))
cat(sprintf("  Columns: %s\n", paste(names(umap_meta), collapse = ", ")))
cat(sprintf("  Cell types: %s\n", paste(unique(umap_meta$cell_type), collapse = ", ")))
cat(sprintf("  Conditions: %s\n", paste(unique(umap_meta$condition), collapse = ", ")))
stopifnot(nrow(umap_meta) > 0)
stopifnot(all(c("UMAP1", "UMAP2", "cell_type", "condition") %in% names(umap_meta)))

# --- Expression matrix ---
cat("\n4. Expression matrix:\n")
cat(sprintf("  Dimensions: %d genes x %d cells\n", nrow(expr_matrix), ncol(expr_matrix)))
cat(sprintf("  Class: %s\n", class(expr_matrix)[1]))
stopifnot(ncol(expr_matrix) == nrow(umap_meta))
stopifnot(inherits(expr_matrix, "dgCMatrix") || inherits(expr_matrix, "matrix"))

# --- Key genes ---
cat("\n5. Key gene check:\n")
key_genes <- c("Myog", "Mstn", "Myh4", "Pax7", "Pecam1", "Pdgfra", "Ptprc")
for (g in key_genes) {
  present <- g %in% gene_list
  in_matrix <- g %in% rownames(expr_matrix)
  cat(sprintf("  %s: gene_list=%s, matrix=%s\n", g, present, in_matrix))
}

# --- DE results ---
cat("\n6. DE results:\n")
cat(sprintf("  Cell types with DE: %s\n", paste(names(de_results), collapse = ", ")))
for (ct in names(de_results)) {
  de <- de_results[[ct]]
  n_sig <- sum(abs(de$avg_log2FC) > 0.25 & de$p_val_adj < 0.05, na.rm = TRUE)
  cat(sprintf("  %s: %d total, %d significant\n", ct, nrow(de), n_sig))
  stopifnot(all(c("gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2") %in% names(de)))
}

# --- Marker genes ---
cat("\n7. Marker genes:\n")
cat(sprintf("  Total markers: %d\n", nrow(marker_genes)))
stopifnot(all(c("gene", "cluster", "avg_log2FC") %in% names(marker_genes)))

# --- Color palette ---
cat("\n8. Color palette:\n")
cat(sprintf("  %d colors defined\n", length(color_palette)))
missing_ct <- setdiff(unique(umap_meta$cell_type), names(color_palette))
if (length(missing_ct) > 0) {
  cat(sprintf("  [WARN] Missing colors for: %s\n", paste(missing_ct, collapse = ", ")))
}

# --- Summary stats ---
cat("\n9. Summary stats:\n")
print(summary_stats)

cat("\n=== Validation complete ===\n")
