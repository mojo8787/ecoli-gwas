#!/usr/bin/env python3
"""
04_build_pa_matrix.py
Build gene presence/absence matrix using one of two strategies:

  FAST PATH  (default, use_bvbrc_pa_matrix: true in config):
    Download pre-computed protein family (plfam_id) assignments from BV-BRC.
    Each plfam_id is a species-level protein family cluster – equivalent to
    a gene cluster produced by Panaroo/Roary but without local compute.

  SLOW PATH  (use_bvbrc_pa_matrix: false):
    Expects Prokka GFF files already in prokka_dir and Panaroo output in
    panaroo_dir (run 04b_run_prokka.py and 05_run_panaroo.py separately).
    Converts Panaroo's gene_presence_absence.csv to a binary matrix.
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml
from loguru import logger
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.yaml"
HEADERS = {"Accept": "application/json"}


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def setup_logging(cfg):
    log_dir = ROOT / cfg["paths"]["logs"]
    log_dir.mkdir(parents=True, exist_ok=True)
    logger.add(log_dir / "04_pa_matrix.log", rotation="10 MB", level="DEBUG")


def load_passing_genomes(cfg: dict) -> list:
    path = ROOT / cfg["paths"]["passing_genomes"]
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run 02_qc_genomes.py first.")
    return path.read_text().strip().splitlines()


# ── Fast path: BV-BRC plfam API ──────────────────────────────────────────────

def fetch_plfam_batch(genome_ids: list, base_url: str) -> pd.DataFrame:
    """Fetch plfam_id assignments for a batch of genome IDs."""
    ids_str = ",".join(genome_ids)
    url = f"{base_url}/genome_feature/"
    params = {
        "in(genome_id,({ids}))".format(ids=ids_str): "",
        "eq(annotation,PATRIC)": "",
        "eq(feature_type,CDS)": "",
        "select": "genome_id,plfam_id",
        "limit": 200000,
    }

    for attempt in range(5):
        try:
            resp = requests.get(url, params=params, headers=HEADERS, timeout=120)
            resp.raise_for_status()
            data = resp.json()
            return pd.DataFrame(data)
        except Exception as e:
            if attempt == 4:
                raise
            wait = 2 ** attempt
            logger.warning(f"Retry {attempt+1}/5: {e}. Waiting {wait}s")
            time.sleep(wait)


def build_pa_matrix_from_bvbrc(genome_ids: list, cfg: dict) -> pd.DataFrame:
    """Download plfam assignments and pivot to presence/absence matrix."""
    base_url = cfg["bvbrc"]["base_url"]
    batch_size = 50  # BV-BRC URL length limit: ~50 genome IDs per request

    all_records = []
    batches = [genome_ids[i:i+batch_size] for i in range(0, len(genome_ids), batch_size)]

    for batch in tqdm(batches, desc="Fetching plfam data"):
        df_batch = fetch_plfam_batch(batch, base_url)
        all_records.append(df_batch)
        time.sleep(0.2)

    df = pd.concat(all_records, ignore_index=True)
    df = df.dropna(subset=["plfam_id"])
    df = df[df["plfam_id"] != ""]
    logger.info(f"Total CDS records: {len(df)}, unique plfam_ids: {df['plfam_id'].nunique()}")

    # Build presence/absence matrix
    pa = pd.crosstab(df["genome_id"], df["plfam_id"])
    pa = pa.clip(upper=1)  # paralogs counted as present once
    logger.info(f"PA matrix shape (before filtering): {pa.shape}")
    return pa


# ── Slow path: Panaroo output ─────────────────────────────────────────────────

def build_pa_matrix_from_panaroo(cfg: dict) -> pd.DataFrame:
    """Parse Panaroo gene_presence_absence.csv into a binary matrix."""
    panaroo_dir = ROOT / cfg["paths"]["panaroo_dir"]
    pa_csv = panaroo_dir / "gene_presence_absence.csv"
    if not pa_csv.exists():
        raise FileNotFoundError(
            f"Panaroo output not found at {pa_csv}. "
            "Run 04b_run_prokka.py and 05_run_panaroo.py first, "
            "or set use_bvbrc_pa_matrix: true in config.yaml."
        )

    df = pd.read_csv(pa_csv, low_memory=False, index_col=0)
    # Columns: 'Non-unique Gene name', 'Annotation', then genome columns
    meta_cols = ["Non-unique Gene name", "Annotation"]
    genome_cols = [c for c in df.columns if c not in meta_cols]
    pa = df[genome_cols].T  # rows=genomes, columns=genes
    pa = (pa != "").astype(int)
    logger.info(f"Panaroo PA matrix shape: {pa.shape}")
    return pa


# ── Shared: filter by allele frequency ────────────────────────────────────────

def filter_by_allele_freq(pa: pd.DataFrame, min_af: float, max_af: float) -> pd.DataFrame:
    """Remove genes that are (near) fixed or (near) absent in the cohort."""
    freq = pa.mean(axis=0)
    before = pa.shape[1]
    pa = pa.loc[:, (freq >= min_af) & (freq <= max_af)]
    after = pa.shape[1]
    logger.info(f"AF filter ({min_af}–{max_af}): kept {after}/{before} gene clusters")
    return pa


def main():
    cfg = load_config()
    setup_logging(cfg)
    logger.info("=== Step 4: Build Gene Presence/Absence Matrix ===")

    genome_ids = load_passing_genomes(cfg)
    logger.info(f"Using {len(genome_ids)} QC-passing genomes")

    use_bvbrc = cfg["cohort"]["use_bvbrc_pa_matrix"]

    if use_bvbrc:
        logger.info("Using FAST PATH: BV-BRC plfam API")
        pa = build_pa_matrix_from_bvbrc(genome_ids, cfg)
    else:
        logger.info("Using SLOW PATH: Panaroo gene_presence_absence.csv")
        pa = build_pa_matrix_from_panaroo(cfg)

    # Ensure only QC-passing genomes are in the matrix
    common = [g for g in genome_ids if g in pa.index]
    pa = pa.loc[common]

    # Filter by allele frequency
    gwas_cfg = cfg["gwas"]
    pa = filter_by_allele_freq(pa, gwas_cfg["min_af"], gwas_cfg["max_af"])

    # Save
    out_path = ROOT / cfg["paths"]["gene_pa_matrix"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pa.to_csv(out_path, sep="\t")
    logger.info(f"Gene PA matrix saved: {pa.shape} → {out_path}")

    # Pangenome stats
    n_genomes = len(pa)
    core_thresh = 0.95
    soft_core_thresh = 0.15
    shell_thresh = 0.05

    freqs = pa.mean(axis=0)
    n_core = (freqs >= core_thresh).sum()
    n_soft = ((freqs >= soft_core_thresh) & (freqs < core_thresh)).sum()
    n_shell = ((freqs >= shell_thresh) & (freqs < soft_core_thresh)).sum()
    n_cloud = (freqs < shell_thresh).sum()

    print("\n=== Pangenome Statistics ===")
    print(f"Genomes          : {n_genomes}")
    print(f"Total gene clusters: {pa.shape[1]}")
    print(f"Core  (>={core_thresh*100:.0f}%)      : {n_core}")
    print(f"Soft-core ({soft_core_thresh*100:.0f}-{core_thresh*100:.0f}%): {n_soft}")
    print(f"Shell ({shell_thresh*100:.0f}-{soft_core_thresh*100:.0f}%)  : {n_shell}")
    print(f"Cloud (<{shell_thresh*100:.0f}%)       : {n_cloud}")

    logger.info("=== Step 4 Complete ===")


if __name__ == "__main__":
    main()
