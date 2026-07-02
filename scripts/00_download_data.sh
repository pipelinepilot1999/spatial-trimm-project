#!/usr/bin/env bash
# 00_download_data.sh — fetch the raw ABC Atlas MERFISH files needed by stages 01 and 04.
# These are gitignored (too large to commit); everything downstream is reproducible once present.
# Source: Allen Brain Cell Atlas public S3 bucket (no credentials required).
# Dataset: MERFISH-C57BL6J-638850 (Zhang et al. 2023), release 20260415 / views 20241115 & 20240330.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw
BASE="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"

echo "[1/2] whole-brain cell metadata + cluster annotation (~1 GB) — used by 01_make_cortical_crop.py"
curl -fL -o data/raw/cell_metadata_with_cluster_annotation.csv \
  "$BASE/metadata/MERFISH-C57BL6J-638850/20241115/views/cell_metadata_with_cluster_annotation.csv"

echo "[2/2] section-.35 expression h5ad (~157 MB) — used by 04_build_de_input.py"
curl -fL -o data/raw/C57BL6J-638850.35-raw.h5ad \
  "$BASE/expression_matrices/MERFISH-C57BL6J-638850-sections/20240330/C57BL6J-638850.35-raw.h5ad"

echo "done. Files in data/raw/:"
ls -lh data/raw/
