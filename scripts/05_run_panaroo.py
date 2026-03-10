#!/usr/bin/env python3
"""
05_run_panaroo.py  (SLOW PATH - optional)
Build pangenome from Prokka GFF files using Panaroo.
Only needed if use_bvbrc_pa_matrix: false in config.yaml.
Requires: panaroo installed (conda install -c bioconda panaroo)
"""

import subprocess
from pathlib import Path

import yaml
from loguru import logger

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.yaml"


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def setup_logging(cfg):
    log_dir = ROOT / cfg["paths"]["logs"]
    log_dir.mkdir(parents=True, exist_ok=True)
    logger.add(log_dir / "05_panaroo.log", rotation="50 MB", level="DEBUG")


def main():
    cfg = load_config()
    setup_logging(cfg)
    logger.info("=== Step 5: Run Panaroo (Slow Path) ===")

    if cfg["cohort"]["use_bvbrc_pa_matrix"]:
        logger.warning(
            "use_bvbrc_pa_matrix is True in config — Panaroo not needed. "
            "Set it to False to use this script. Exiting."
        )
        return

    prokka_dir = ROOT / cfg["paths"]["prokka_dir"]
    panaroo_dir = ROOT / cfg["paths"]["panaroo_dir"]
    panaroo_dir.mkdir(parents=True, exist_ok=True)

    # Collect all GFF files
    gff_files = sorted(prokka_dir.glob("*/*.gff"))
    if not gff_files:
        logger.error(f"No GFF files found in {prokka_dir}. Run 04b_run_prokka.py first.")
        return

    logger.info(f"Found {len(gff_files)} GFF files")

    # Check if already done
    pa_csv = panaroo_dir / "gene_presence_absence.csv"
    if pa_csv.exists():
        logger.info(f"Panaroo output already exists at {pa_csv}. Delete to re-run.")
        return

    threads = cfg["compute"]["panaroo_threads"]
    cmd = [
        "panaroo",
        "-i", *[str(g) for g in gff_files],
        "-o", str(panaroo_dir),
        "--clean-mode", "strict",
        "-t", str(threads),
        "--remove-invalid-genes",
    ]

    logger.info(f"Running Panaroo with {threads} threads...")
    logger.info(f"Command: {' '.join(cmd[:5])} ... [{len(gff_files)} GFF files]")

    log_file = ROOT / cfg["paths"]["logs"] / "panaroo_stdout.log"
    with open(log_file, "w") as lf:
        result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        logger.error(f"Panaroo failed (exit {result.returncode}). Check {log_file}")
        return

    logger.info(f"Panaroo complete. Output at {panaroo_dir}")
    logger.info("=== Step 5 Complete ===")


if __name__ == "__main__":
    main()
