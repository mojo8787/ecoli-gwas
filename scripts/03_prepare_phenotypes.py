#!/usr/bin/env python3
"""
03_prepare_phenotypes.py
Build a clean binary phenotype matrix (genomes x antibiotics).
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from loguru import logger

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.yaml"

# Antibiotic name synonyms → canonical name
SYNONYMS = {
    "amoxicillin": "ampicillin",
    "amoxicillin/clavulanic acid": "amoxicillin-clavulanate",
    "amoxicillin-clavulanate": "amoxicillin-clavulanate",
    "trimethoprim-sulfamethoxazole": "trimethoprim",
    "co-trimoxazole": "trimethoprim",
    "sulfamethoxazole/trimethoprim": "trimethoprim",
    "ceftriaxone": "ceftriaxone",
    "cefotaxime": "ceftriaxone",
    "imipenem": "meropenem",
    "ertapenem": "meropenem",
    "colistin": "polymyxin b",
}


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def setup_logging(cfg):
    log_dir = ROOT / cfg["paths"]["logs"]
    log_dir.mkdir(parents=True, exist_ok=True)
    logger.add(log_dir / "03_phenotypes.log", rotation="10 MB", level="DEBUG")


def normalize_antibiotic(name: str) -> str:
    n = name.lower().strip()
    return SYNONYMS.get(n, n)


def load_amr_data(cfg: dict) -> pd.DataFrame:
    amr_path = ROOT / cfg["paths"]["raw_amr"] / "amr_phenotypes.json"
    if not amr_path.exists():
        raise FileNotFoundError(f"AMR data not found: {amr_path}. Run 01_download_data.py first.")
    with open(amr_path) as f:
        data = json.load(f)
    return pd.DataFrame(data)


def load_passing_genomes(cfg: dict) -> set:
    path = ROOT / cfg["paths"]["passing_genomes"]
    if not path.exists():
        raise FileNotFoundError(f"Passing genomes not found: {path}. Run 02_qc_genomes.py first.")
    return set(path.read_text().strip().splitlines())


def resolve_conflicts(group: pd.DataFrame) -> float:
    """
    Given multiple AMR calls for one genome+antibiotic, return 1, 0, or NaN.
    Uses majority vote; returns NaN if tied.
    """
    r_count = (group["resistant_phenotype"] == "Resistant").sum()
    s_count = (group["resistant_phenotype"] == "Susceptible").sum()
    if r_count > s_count:
        return 1.0
    elif s_count > r_count:
        return 0.0
    return np.nan


def build_phenotype_matrix(amr_df: pd.DataFrame, passing_ids: set, cfg: dict) -> pd.DataFrame:
    """Create a binary (genome x antibiotic) matrix."""
    ph = cfg["phenotypes"]

    # Filter to QC-passing genomes
    amr_df = amr_df[amr_df["genome_id"].isin(passing_ids)].copy()
    logger.info(f"AMR records after QC filter: {len(amr_df)}")

    # Normalize antibiotic names
    amr_df["antibiotic"] = amr_df["antibiotic"].apply(normalize_antibiotic)

    # Filter to antibiotics of interest
    ab_interest = [a.lower().strip() for a in ph["antibiotics_of_interest"]]
    if ab_interest:
        amr_df = amr_df[amr_df["antibiotic"].isin(ab_interest)]
    logger.info(f"AMR records after antibiotic filter: {len(amr_df)}")

    # Resolve conflicts (majority vote per genome+antibiotic)
    resolved = (
        amr_df.groupby(["genome_id", "antibiotic"])
        .apply(resolve_conflicts)
        .reset_index()
    )
    resolved.columns = ["genome_id", "antibiotic", "phenotype"]

    # Pivot to wide format
    matrix = resolved.pivot(index="genome_id", columns="antibiotic", values="phenotype")

    # Filter antibiotics with enough cases
    min_r = ph["min_resistant"]
    min_s = ph["min_susceptible"]
    eligible_ab = []
    for ab in matrix.columns:
        col = matrix[ab].dropna()
        n_r = (col == 1).sum()
        n_s = (col == 0).sum()
        if n_r >= min_r and n_s >= min_s:
            eligible_ab.append(ab)
        else:
            logger.info(f"  Dropping {ab}: R={n_r}, S={n_s} (need >={min_r} R and >={min_s} S)")

    matrix = matrix[eligible_ab]
    logger.info(f"Eligible antibiotics: {eligible_ab}")
    return matrix


def main():
    cfg = load_config()
    setup_logging(cfg)
    logger.info("=== Step 3: Prepare Phenotypes ===")

    amr_df = load_amr_data(cfg)
    passing_ids = load_passing_genomes(cfg)
    logger.info(f"Passing genomes: {len(passing_ids)}")

    matrix = build_phenotype_matrix(amr_df, passing_ids, cfg)

    out_path = ROOT / cfg["paths"]["phenotype_matrix"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(out_path, sep="\t")
    logger.info(f"Phenotype matrix saved: {matrix.shape} → {out_path}")

    # Summary
    print("\n=== Phenotype Summary ===")
    print(f"Genomes in matrix : {len(matrix)}")
    print(f"Antibiotics       : {list(matrix.columns)}")
    print()
    print(f"{'Antibiotic':<30} {'Resistant':>10} {'Susceptible':>12} {'Missing':>8} {'Balance':>8}")
    print("-" * 72)
    for ab in matrix.columns:
        col = matrix[ab]
        n_r = int((col == 1).sum())
        n_s = int((col == 0).sum())
        n_na = int(col.isna().sum())
        balance = n_r / (n_r + n_s) if (n_r + n_s) > 0 else float("nan")
        print(f"{ab:<30} {n_r:>10} {n_s:>12} {n_na:>8} {balance:>8.2f}")

    logger.info("=== Step 3 Complete ===")


if __name__ == "__main__":
    main()
