#!/usr/bin/env python3
"""
02_qc_genomes.py
Compute assembly statistics with Biopython and filter poor-quality genomes.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio import SeqIO
from loguru import logger
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.yaml"


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def setup_logging(cfg):
    log_dir = ROOT / cfg["paths"]["logs"]
    log_dir.mkdir(parents=True, exist_ok=True)
    logger.add(log_dir / "02_qc.log", rotation="10 MB", level="DEBUG")


def compute_n50(lengths: list[int]) -> tuple[int, int]:
    """Return (N50, L50) for a list of contig lengths."""
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    cumsum = 0
    for i, length in enumerate(sorted_lengths):
        cumsum += length
        if cumsum >= total / 2:
            return length, i + 1
    return 0, 0


def assess_genome(fasta_path: Path) -> dict:
    """Compute QC metrics for a single genome FASTA."""
    lengths = []
    gc_counts = []
    n_counts = []

    try:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            seq = str(rec.seq).upper()
            lengths.append(len(seq))
            gc = seq.count("G") + seq.count("C")
            gc_counts.append(gc)
            n_counts.append(seq.count("N"))
    except TimeoutError as e:
        logger.warning(f"Timeout while reading {fasta_path}: {e}")
        return {"genome_id": fasta_path.stem, "error": "read_timeout"}
    except OSError as e:
        logger.warning(f"I/O error while reading {fasta_path}: {e}")
        return {"genome_id": fasta_path.stem, "error": "read_error"}

    if not lengths:
        return {"genome_id": fasta_path.stem, "error": "empty_file"}

    total_length = sum(lengths)
    total_gc = sum(gc_counts)
    total_n = sum(n_counts)
    n50, l50 = compute_n50(lengths)
    contigs_over_500 = sum(1 for l in lengths if l >= 500)

    return {
        "genome_id": fasta_path.stem,
        "total_length": total_length,
        "num_contigs": len(lengths),
        "contigs_ge_500bp": contigs_over_500,
        "n50": n50,
        "l50": l50,
        "gc_content": total_gc / max(total_length, 1),
        "n_bases": total_n,
        "n_fraction": total_n / max(total_length, 1),
        "max_contig": max(lengths),
        "min_contig": min(lengths),
    }


def apply_filters(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Add PASS/FAIL column and fail_reason to the QC dataframe."""
    qc = cfg["qc"]
    reasons = []

    for _, row in df.iterrows():
        r = []
        if row.get("error"):
            r.append(row["error"])
        else:
            if row["total_length"] < qc["min_length"]:
                r.append(f"length_too_short({row['total_length']})")
            if row["total_length"] > qc["max_length"]:
                r.append(f"length_too_long({row['total_length']})")
            if row["num_contigs"] > qc["max_contigs"]:
                r.append(f"too_many_contigs({row['num_contigs']})")
            if row["n50"] < qc["min_n50"]:
                r.append(f"low_n50({row['n50']})")
            if row["gc_content"] < qc["min_gc"]:
                r.append(f"low_gc({row['gc_content']:.3f})")
            if row["gc_content"] > qc["max_gc"]:
                r.append(f"high_gc({row['gc_content']:.3f})")
        # Ensure all reasons are strings and skip nulls
        reasons.append("; ".join(str(reason) for reason in r if pd.notna(reason)) if r else "")

    df = df.copy()
    df["fail_reason"] = reasons
    df["qc_pass"] = df["fail_reason"] == ""
    return df


def main():
    cfg = load_config()
    setup_logging(cfg)
    logger.info("=== Step 2: QC Genomes ===")

    genome_dir = ROOT / cfg["paths"]["raw_genomes"]
    fasta_files = sorted(genome_dir.glob("*.fna"))
    logger.info(f"Found {len(fasta_files)} FASTA files to assess")

    if not fasta_files:
        logger.error(f"No .fna files found in {genome_dir}. Run 01_download_data.py first.")
        return

    # Compute stats
    records = []
    for fasta in tqdm(fasta_files, desc="Computing QC stats"):
        records.append(assess_genome(fasta))

    df = pd.DataFrame(records)
    df = apply_filters(df, cfg)

    # Save full report
    qc_report_path = ROOT / cfg["paths"]["qc_report"]
    qc_report_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(qc_report_path, sep="\t", index=False)
    logger.info(f"QC report saved to {qc_report_path}")

    # Save passing genome IDs
    passing = df[df["qc_pass"]]["genome_id"].tolist()
    passing_path = ROOT / cfg["paths"]["passing_genomes"]
    passing_path.parent.mkdir(parents=True, exist_ok=True)
    passing_path.write_text("\n".join(passing) + "\n")
    logger.info(f"Passing genomes: {len(passing)} / {len(df)}")

    # Summary by fail reason
    failed = df[~df["qc_pass"]]
    if len(failed) > 0:
        logger.info("Failure breakdown:")
        for reason_group in failed["fail_reason"].str.split(";").explode().str.strip().value_counts().items():
            logger.info(f"  {reason_group[0]}: {reason_group[1]}")

    # Print summary table
    print("\n=== QC Summary ===")
    print(f"Total genomes assessed : {len(df)}")
    print(f"Passed QC              : {len(passing)}")
    print(f"Failed QC              : {len(failed)}")
    print(f"\nPassing genome stats:")
    passing_df = df[df["qc_pass"]]
    if not passing_df.empty:
        for col in ["total_length", "num_contigs", "n50", "gc_content"]:
            print(f"  {col:20s} mean={passing_df[col].mean():.1f}  min={passing_df[col].min():.1f}  max={passing_df[col].max():.1f}")

    logger.info("=== Step 2 Complete ===")


if __name__ == "__main__":
    main()
