#!/usr/bin/env bash
# Download processed h5 files from GEO GSE272085
# These are Cell Ranger ARC outputs (multiome: GEX + ATAC)
# We extract only the Gene Expression modality in the analysis pipeline

set -euo pipefail

OUTDIR="$(cd "$(dirname "$0")" && pwd)"

echo "Downloading Control (WT) h5 file..."
curl -L -o "${OUTDIR}/GSM8392727_WT_filtered_feature_bc_matrix.h5" \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8392727&format=file&file=GSM8392727%5FWT%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"

echo "Downloading KIC (Cachexia) h5 file..."
curl -L -o "${OUTDIR}/GSM8392728_Cachexia_filtered_feature_bc_matrix.h5" \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8392728&format=file&file=GSM8392728%5FCachexia%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"

echo "Download complete."
ls -lh "${OUTDIR}"/*.h5
