# 01_qc_and_preprocessing.R
# Read multiome h5 files, extract Gene Expression, QC filter, normalize

library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)

set.seed(42)

raw_dir <- "data/raw"
out_dir <- "data/processed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("results/qc", showWarnings = FALSE, recursive = TRUE)

# --- Read h5 files (multiome: extract Gene Expression only) ---
cat("Reading Control h5...\n")
wt_data <- Read10X_h5(
  file.path(raw_dir, "GSM8392727_WT_filtered_feature_bc_matrix.h5")
)
# Multiome h5 returns a list; extract Gene Expression
if (is.list(wt_data)) {
  wt_counts <- wt_data[["Gene Expression"]]
} else {
  wt_counts <- wt_data
}

cat("Reading KIC h5...\n")
kic_data <- Read10X_h5(
  file.path(raw_dir, "GSM8392728_Cachexia_filtered_feature_bc_matrix.h5")
)
if (is.list(kic_data)) {
  kic_counts <- kic_data[["Gene Expression"]]
} else {
  kic_counts <- kic_data
}

cat(sprintf("Control: %d genes x %d nuclei\n", nrow(wt_counts), ncol(wt_counts)))
cat(sprintf("KIC: %d genes x %d nuclei\n", nrow(kic_counts), ncol(kic_counts)))

# --- Create Seurat objects ---
wt <- CreateSeuratObject(counts = wt_counts, project = "Control",
                         min.cells = 3, min.features = 200)
wt$condition <- "Control"

kic <- CreateSeuratObject(counts = kic_counts, project = "KIC",
                          min.cells = 3, min.features = 200)
kic$condition <- "KIC"

# --- QC metrics ---
wt[["percent.mt"]] <- PercentageFeatureSet(wt, pattern = "^mt-")
kic[["percent.mt"]] <- PercentageFeatureSet(kic, pattern = "^mt-")

# QC plots before filtering
p_wt <- VlnPlot(wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3, pt.size = 0) + plot_annotation(title = "Control - Pre-filter")
p_kic <- VlnPlot(kic, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, pt.size = 0) + plot_annotation(title = "KIC - Pre-filter")

ggsave("results/qc/qc_prefilter_control.pdf", p_wt, width = 12, height = 5)
ggsave("results/qc/qc_prefilter_kic.pdf", p_kic, width = 12, height = 5)

# --- Filter ---
cat("Filtering nuclei...\n")
wt <- subset(wt, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
               nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 20)
kic <- subset(kic, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 20)

cat(sprintf("After filtering - Control: %d nuclei\n", ncol(wt)))
cat(sprintf("After filtering - KIC: %d nuclei\n", ncol(kic)))

# QC plots after filtering
p_wt2 <- VlnPlot(wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, pt.size = 0) + plot_annotation(title = "Control - Post-filter")
p_kic2 <- VlnPlot(kic, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0) + plot_annotation(title = "KIC - Post-filter")

ggsave("results/qc/qc_postfilter_control.pdf", p_wt2, width = 12, height = 5)
ggsave("results/qc/qc_postfilter_kic.pdf", p_kic2, width = 12, height = 5)

# --- Normalize with SCTransform v2 ---
cat("Running SCTransform on Control...\n")
wt <- SCTransform(wt, vst.flavor = "v2", verbose = FALSE)

cat("Running SCTransform on KIC...\n")
kic <- SCTransform(kic, vst.flavor = "v2", verbose = FALSE)

# --- Save ---
saveRDS(wt, file.path(out_dir, "wt_sct.rds"))
saveRDS(kic, file.path(out_dir, "kic_sct.rds"))

cat("QC and preprocessing complete.\n")
