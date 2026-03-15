# 04_annotation.R
# Cell type annotation — revised with denervated/catabolic myonuclei
# Based on Zhang et al. 2024 (PMID 39116208) cell type definitions

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(42)

out_dir <- "data/processed"
dir.create("results/annotation", showWarnings = FALSE, recursive = TRUE)

# --- Load ---
cat("Loading clustered object...\n")
merged <- readRDS(file.path(out_dir, "merged_clustered.rds"))
Idents(merged) <- "SCT_snn_res.0.5"

# --- Expanded marker gene panels ---
# Includes denervated/catabolic markers from Zhang et al.
markers <- list(
  "Type IIB myonuclei" = c("Myh4", "Mybpc2", "Tnnt3", "Ttn", "Mylpf"),
  "Type IIX myonuclei" = c("Myh1", "Casq1", "Atp2a1"),
  "Type IIA myonuclei" = c("Myh2", "Tnni1", "Ppara"),
  "Type I myonuclei"   = c("Myh7", "Tnnc1", "Myl3", "Atp2a2"),
  "Denervated myonuclei" = c("Ncam1", "Runx1", "Gadd45a", "Myog", "Chrna1"),
  "Catabolic myonuclei"  = c("Fbxo32", "Trim63"),
  "Satellite cells"    = c("Pax7", "Myf5", "Calcr", "Cd34", "Cxcr4"),
  "Endothelial"        = c("Pecam1", "Cdh5", "Flt1", "Kdr", "Emcn"),
  "Lymphatic endothelial" = c("Flt4", "Mmrn1", "Ccl21a", "Reln"),
  "Fibroblasts/FAPs"   = c("Pdgfra", "Dcn", "Col1a1", "Col3a1"),
  "Tenocytes"          = c("Mkx", "Fmod", "Col11a1", "Tnmd"),
  "Immune"             = c("Ptprc", "Cd68", "Adgre1", "Cd163", "Mrc1"),
  "Pericytes"          = c("Rgs5", "Acta2", "Pdgfrb", "Myh11", "Notch3"),
  "Adipocytes"         = c("Adipoq", "Pparg", "Lep", "Fabp4", "Plin1"),
  "NMJ/MTJ nuclei"    = c("Col22a1", "Chrne", "Musk")
)

all_markers <- unique(unlist(markers))
available_markers <- all_markers[all_markers %in% rownames(merged)]
cat(sprintf("Using %d of %d marker genes present in dataset\n",
            length(available_markers), length(all_markers)))

# --- FindAllMarkers for annotation support ---
cat("Running FindAllMarkers...\n")
merged <- PrepSCTFindMarkers(merged)
all_de <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25,
                         logfc.threshold = 0.25, verbose = FALSE)

top_markers_per_cluster <- all_de %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top_markers_per_cluster, "results/annotation/top_markers_per_cluster.csv",
          row.names = FALSE)

# --- Manual annotation based on marker expression + condition breakdown ---
# Determined by examining:
#   1) Top DE markers per cluster
#   2) Expression of canonical + denervation/catabolic markers
#   3) Condition composition (ctrl vs KIC enrichment)
#   4) Concordance with Zhang et al. cell type definitions

cluster_annotations <- c(
  "0"  = "Fibroblasts/FAPs",        # Pdgfra+, Dcn+, Col1a1+, Col3a1+
  "1"  = "Type IIB myonuclei",      # Myh4 99%, Mybpc2 86%, Ckm 92%, ctrl-dominant
  "2"  = "Denervated myonuclei",    # Ncam1 92%, Runx1 93%, Gadd45a 73%, Myog 25%, KIC-only
  "3"  = "Catabolic myonuclei",     # Fbxo32 100% (avg 4.0), Trim63 100% (avg 4.14), low Myh4, KIC-only
  "4"  = "Type IIB myonuclei",      # Myh4 87%, Ckm 80%, KIC-dominant (stressed IIB)
  "5"  = "Type IIX myonuclei",      # Myh1 56%, Ankrd2+, KIC-dominant
  "6"  = "Type IIA myonuclei",      # Myh2 54%, Myh1 85%, Ppara+, ctrl-dominant
  "7"  = "Pericytes",               # Rgs5+, Myh11+, Gucy1a1+
  "8"  = "NMJ/MTJ nuclei",          # Col22a1 top marker, Ncam1 41%, small mixed cluster
  "9"  = "Immune",                   # Ptprc+, Cd68+, F4/80+
  "10" = "Type IIB myonuclei",      # Mylpf 95%, Myh4 100%, Ckm 99%, KIC-only
  "11" = "Type I myonuclei",        # Myh7 90%, Tnnc1+, Myl3+
  "12" = "Tenocytes",               # Col11a1+, Fmod+, Mkx+ (tendon TF)
  "13" = "Endothelial",             # Pecam1+, Sox17+, Ptprb+
  "14" = "Satellite cells",         # Calcr+, Fgfr4+, Pax7+
  "15" = "Lymphatic endothelial",   # Mmrn1+, Ccl21a+, Flt4+, Reln+
  "16" = "Adipocytes"               # Pck1+, Cfd+, Cidec+, Adipoq+
)

cat("\nCluster annotations:\n")
for (cl in names(cluster_annotations)) {
  cat(sprintf("  Cluster %s -> %s\n", cl, cluster_annotations[cl]))
}

# Apply annotations
cell_type_vec <- unname(cluster_annotations[as.character(Idents(merged))])
merged <- AddMetaData(merged, metadata = cell_type_vec, col.name = "cell_type")

# --- DotPlot of all markers by cell type ---
Idents(merged) <- "cell_type"
p_dot_ct <- DotPlot(merged, features = available_markers, cluster.idents = TRUE) +
  RotatedAxis() +
  ggtitle("Marker Gene Expression by Cell Type") +
  theme(axis.text.x = element_text(size = 8))

ggsave("results/annotation/dotplot_markers_annotated.pdf", p_dot_ct, width = 22, height = 10)

# --- DotPlot by cluster (for reference) ---
Idents(merged) <- "SCT_snn_res.0.5"
p_dot_cl <- DotPlot(merged, features = available_markers, cluster.idents = TRUE) +
  RotatedAxis() +
  ggtitle("Marker Gene Expression by Cluster") +
  theme(axis.text.x = element_text(size = 8))

ggsave("results/annotation/dotplot_markers_by_cluster.pdf", p_dot_cl, width = 22, height = 10)

# --- Annotated UMAPs ---
Idents(merged) <- "cell_type"

p_annotated <- DimPlot(merged, group.by = "cell_type", label = TRUE, repel = TRUE,
                       label.size = 3.5) +
  ggtitle("Cell Type Annotations") + NoLegend()
p_condition <- DimPlot(merged, group.by = "condition") +
  ggtitle("By Condition")

ggsave("results/annotation/umap_annotated.pdf", p_annotated + p_condition,
       width = 16, height = 7)

# --- FeaturePlots for key markers including denervation/catabolic ---
key_markers <- c("Myh4", "Myh1", "Myh2", "Myh7",     # fiber types
                 "Ncam1", "Runx1", "Gadd45a", "Myog",  # denervation
                 "Fbxo32", "Trim63",                    # catabolic/atrophy
                 "Pax7", "Pecam1", "Pdgfra", "Ptprc",  # non-muscle
                 "Rgs5", "Adipoq", "Col22a1")           # other

key_available <- key_markers[key_markers %in% rownames(merged)]
p_feat <- FeaturePlot(merged, features = key_available, ncol = 4,
                      order = TRUE, cols = c("lightgrey", "darkred"))

ggsave("results/annotation/featureplot_key_markers.pdf", p_feat,
       width = 20, height = 20)

# --- Annotation summary ---
summary_table <- merged@meta.data %>%
  group_by(cell_type, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0)

cat("\nCell type composition:\n")
print(as.data.frame(summary_table))
write.csv(summary_table, "results/annotation/cell_type_summary.csv", row.names = FALSE)

# --- Heatmap of marker expression ---
top5 <- all_de %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  pull(gene) %>%
  unique()

# Use ScaleData for heatmap
Idents(merged) <- "cell_type"
p_heat <- DoHeatmap(merged, features = top5, group.by = "cell_type",
                    size = 3) + NoLegend()
ggsave("results/annotation/heatmap_top_markers.pdf", p_heat, width = 16, height = 12)

# --- Save ---
saveRDS(merged, file.path(out_dir, "merged_annotated.rds"))
saveRDS(all_de, file.path(out_dir, "all_cluster_markers.rds"))

cat("Annotation complete.\n")
