# 05_differential_expression.R
# Per cell-type DE: KIC vs Control
# Uses RNA assay with log-normalization (avoids SCT multi-model issues with subset)

library(Seurat)
library(dplyr)

set.seed(42)

out_dir <- "data/processed"
dir.create("results/de", showWarnings = FALSE, recursive = TRUE)

# --- Load ---
cat("Loading annotated object...\n")
merged <- readRDS(file.path(out_dir, "merged_annotated.rds"))

# Switch to RNA assay and normalize for DE
# (Seurat recommends RNA assay for cross-condition DE to avoid SCT multi-model issues)
DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged, verbose = FALSE)

# Set active identity to cell_type
Idents(merged) <- "cell_type"

# --- Per cell-type DE ---
cell_types <- unique(merged$cell_type)
cell_types <- cell_types[!is.na(cell_types) & cell_types != "Unknown"]

de_results <- list()

for (ct in cell_types) {
  cat(sprintf("DE for %s...\n", ct))

  cells_ct <- WhichCells(merged, idents = ct)

  # Need cells from both conditions
  conditions_present <- unique(merged$condition[cells_ct])
  if (length(conditions_present) < 2) {
    cat(sprintf("  Skipping %s: only one condition present\n", ct))
    next
  }

  n_kic <- sum(merged$condition[cells_ct] == "KIC")
  n_ctrl <- sum(merged$condition[cells_ct] == "Control")
  cat(sprintf("  KIC: %d, Control: %d nuclei\n", n_kic, n_ctrl))

  if (n_kic < 10 || n_ctrl < 10) {
    cat(sprintf("  Skipping %s: too few cells in one condition\n", ct))
    next
  }

  tryCatch({
    # Subset to this cell type, then run DE by condition
    sub_obj <- subset(merged, idents = ct)
    Idents(sub_obj) <- "condition"

    de <- FindMarkers(
      sub_obj,
      ident.1 = "KIC",
      ident.2 = "Control",
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.1,
      verbose = FALSE
    )

    de$gene <- rownames(de)
    de$cell_type <- ct
    de_results[[ct]] <- de

    # Save individual CSV
    ct_clean <- gsub("[/ ]", "_", ct)
    write.csv(de, sprintf("results/de/de_%s.csv", ct_clean), row.names = FALSE)

    n_sig <- sum(abs(de$avg_log2FC) > 0.25 & de$p_val_adj < 0.05, na.rm = TRUE)
    cat(sprintf("  %d significant DE genes (|log2FC| > 0.25, padj < 0.05)\n", n_sig))

  }, error = function(e) {
    cat(sprintf("  Error in %s: %s\n", ct, e$message))
  })
}

# --- Validate key genes from paper ---
key_genes <- c("Myog", "Mstn", "Chrna1", "Chrne", "Fbxo32", "Trim63")
cat("\n=== Key gene validation ===\n")

for (ct in names(de_results)) {
  de <- de_results[[ct]]
  found <- de %>% filter(gene %in% key_genes & p_val_adj < 0.05)
  if (nrow(found) > 0) {
    cat(sprintf("\n%s:\n", ct))
    for (i in seq_len(nrow(found))) {
      cat(sprintf("  %s: log2FC = %.3f, padj = %.2e\n",
                  found$gene[i], found$avg_log2FC[i], found$p_val_adj[i]))
    }
  }
}

# --- Save combined results ---
saveRDS(de_results, file.path(out_dir, "de_results.rds"))

# Combined CSV
all_de <- bind_rows(de_results)
write.csv(all_de, "results/de/de_all_celltypes.csv", row.names = FALSE)

cat("\nDifferential expression complete.\n")
