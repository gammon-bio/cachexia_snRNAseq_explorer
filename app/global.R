# global.R
# Load pre-computed data at app startup (no Seurat dependency)

library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(scales)
library(waiter)
library(htmltools)

# --- Load pre-computed data ---
data_dir <- file.path("data")

umap_meta    <- readRDS(file.path(data_dir, "umap_metadata.rds"))
expr_matrix  <- readRDS(file.path(data_dir, "expression_matrix.rds"))
de_results   <- readRDS(file.path(data_dir, "de_results.rds"))
marker_genes <- readRDS(file.path(data_dir, "marker_genes.rds"))
gene_list    <- readRDS(file.path(data_dir, "gene_list.rds"))
color_palette <- readRDS(file.path(data_dir, "color_palette.rds"))
summary_stats <- readRDS(file.path(data_dir, "summary_stats.rds"))

# Derived values
cell_types <- sort(unique(umap_meta$cell_type))
conditions <- sort(unique(umap_meta$condition))
n_total <- nrow(umap_meta)

# Source modules
source("modules/mod_gene_explorer.R")
source("modules/mod_de_results.R")
source("modules/mod_cell_composition.R")
source("modules/mod_marker_genes.R")
