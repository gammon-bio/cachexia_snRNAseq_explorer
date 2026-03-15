# plot_annotated_dotplot.R
# DotPlot of canonical marker genes with cell type annotations on the y-axis

library(Seurat)
library(ggplot2)

set.seed(42)

dir.create("results/annotation", showWarnings = FALSE, recursive = TRUE)

# --- Load annotated object ---
cat("Loading annotated Seurat object...\n")
merged <- readRDS("data/processed/merged_annotated.rds")

# Set identities to cell type names
Idents(merged) <- "cell_type"

# --- Marker gene panels (same as 04_annotation.R) ---
markers <- list(
  "Type IIB myonuclei" = c("Myh4", "Mybpc2", "Tnnt3", "Ttn", "Mylpf"),
  "Type IIX myonuclei" = c("Myh1", "Casq1", "Atp2a1"),
  "Type IIA myonuclei" = c("Myh2", "Tnni1"),
  "Type I myonuclei"   = c("Myh7", "Tnnc1", "Myl3", "Atp2a2"),
  "Satellite cells"    = c("Pax7", "Myf5", "Calcr", "Cd34", "Cxcr4"),
  "Endothelial"        = c("Pecam1", "Cdh5", "Flt1", "Kdr", "Emcn"),
  "Fibroblasts/FAPs"   = c("Pdgfra", "Dcn", "Col1a1", "Col3a1"),
  "Immune"             = c("Ptprc", "Cd68", "Adgre1", "Cd163", "Mrc1"),
  "Pericytes"          = c("Rgs5", "Acta2", "Pdgfrb", "Myh11", "Notch3"),
  "Adipocytes"         = c("Adipoq", "Pparg", "Lep", "Fabp4", "Plin1"),
  "NMJ/MTJ"            = c("Chrna1", "Chrne", "Musk", "Tnmd", "Col22a1")
)

all_markers <- unique(unlist(markers))

# Filter to genes present in the dataset
available_markers <- all_markers[all_markers %in% rownames(merged)]
cat(sprintf("Using %d of %d marker genes present in dataset\n",
            length(available_markers), length(all_markers)))

# --- DotPlot ---
cat("Creating DotPlot...\n")
p_dot <- DotPlot(merged, features = available_markers, cluster.idents = TRUE) +
  RotatedAxis() +
  ggtitle("Marker Gene Expression by Cell Type") +
  theme(axis.text.x = element_text(size = 8))

out_file <- "results/annotation/dotplot_markers_annotated.pdf"
ggsave(out_file, p_dot, width = 20, height = 8)
cat(sprintf("Saved to %s\n", out_file))

cat("Done.\n")
