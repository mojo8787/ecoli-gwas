#!/usr/bin/env bash
# ============================================================
# run_pipeline.sh — Full E. coli GWAS pipeline
# Usage: bash run_pipeline.sh
# ============================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "$0")/scripts" && pwd)"
PYTHON="$(cd "$(dirname "$0")" && pwd)/venv/bin/python"

echo "============================================"
echo "  E. coli GWAS Pipeline"
echo "============================================"

echo ""
echo "[1/6] Downloading genomes + AMR phenotypes..."
$PYTHON "$SCRIPTS_DIR/01_download_data.py"

echo ""
echo "[2/6] QC genomes..."
$PYTHON "$SCRIPTS_DIR/02_qc_genomes.py"

echo ""
echo "[3/6] Preparing phenotype matrix..."
$PYTHON "$SCRIPTS_DIR/03_prepare_phenotypes.py"

echo ""
echo "[4/6] Building gene presence/absence matrix..."
$PYTHON "$SCRIPTS_DIR/04_build_pa_matrix.py"

echo ""
echo "[5/6] Running GWAS..."
$PYTHON "$SCRIPTS_DIR/06_run_gwas.py"

echo ""
echo "[6/6] Generating visualizations..."
$PYTHON "$SCRIPTS_DIR/07_visualize.py"

echo ""
echo "============================================"
echo "  Pipeline complete!"
echo "  Results: results/figures/"
echo "  Hits:    results/significant_hits_summary.tsv"
echo "============================================"
