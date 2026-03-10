#!/usr/bin/env python3
"""
01_download_data.py
Download E. coli genomes and AMR phenotypes from BV-BRC public API.
"""

import json
import time
from pathlib import Path

import pandas as pd
import requests
import yaml
from loguru import logger
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.yaml"

SESSION = requests.Session()
SESSION.headers.update({"Accept": "application/json"})


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def setup_logging(cfg):
    log_dir = ROOT / cfg["paths"]["logs"]
    log_dir.mkdir(parents=True, exist_ok=True)
    logger.add(log_dir / "01_download.log", rotation="10 MB", level="DEBUG")


# ── BV-BRC RQL helpers ────────────────────────────────────────────────────────

def bvbrc_query(endpoint: str, rql_conditions: list, select_fields: list,
                limit: int = 5000, retries: int = 5, timeout: int = 90) -> list:
    """
    Query BV-BRC API using proper RQL syntax.
    Builds the query string manually to avoid URL-encoding the parentheses.
    """
    base = "https://www.bv-brc.org/api"
    url = f"{base}/{endpoint}/"
    all_records = []
    offset = 0

    while True:
        rql_parts = rql_conditions + [
            f"select({','.join(select_fields)})",
            f"limit({limit},{offset})",
        ]
        query_string = "&".join(rql_parts)
        full_url = f"{url}?{query_string}"

        for attempt in range(retries):
            try:
                resp = SESSION.get(full_url, timeout=timeout)
                resp.raise_for_status()
                data = resp.json()
                break
            except (requests.exceptions.JSONDecodeError, requests.exceptions.Timeout) as e:
                if attempt == retries - 1:
                    raise
                wait = 2 ** attempt
                logger.warning(f"Retry {attempt+1}/{retries}: {e}. Waiting {wait}s")
                time.sleep(wait)
            except Exception as e:
                if attempt == retries - 1:
                    raise
                wait = 2 ** attempt
                logger.warning(f"Retry {attempt+1}/{retries}: {e}. Waiting {wait}s")
                time.sleep(wait)

        if not data:
            break
        all_records.extend(data)
        logger.debug(f"  offset={offset} fetched={len(data)} total={len(all_records)}")

        content_range = resp.headers.get("Content-Range", "")
        if content_range:
            try:
                total = int(content_range.split("/")[-1])
                if offset + limit >= total:
                    break
            except ValueError:
                pass

        if len(data) < limit:
            break
        offset += limit
        time.sleep(0.2)

    return all_records


def bvbrc_query_resumable(endpoint: str, rql_conditions: list, select_fields: list,
                          limit: int, partial_path: Path, retries: int = 5,
                          timeout: int = 120, save_every_batches: int = 20) -> list:
    """
    Like bvbrc_query but saves progress to partial_path so we can resume after timeout/JSON errors.
    """
    base = "https://www.bv-brc.org/api"
    url = f"{base}/{endpoint}/"

    if partial_path.exists():
        try:
            with open(partial_path) as f:
                all_records = json.load(f)
            offset = len(all_records)
            logger.info(f"Resuming AMR fetch from offset={offset} ({len(all_records)} records already saved)")
        except (json.JSONDecodeError, Exception) as e:
            logger.warning(f"Could not load partial file: {e}. Starting from 0.")
            all_records = []
            offset = 0
    else:
        all_records = []
        offset = 0

    batch_count = 0

    while True:
        rql_parts = rql_conditions + [
            f"select({','.join(select_fields)})",
            f"limit({limit},{offset})",
        ]
        query_string = "&".join(rql_parts)
        full_url = f"{url}?{query_string}"

        for attempt in range(retries):
            try:
                resp = SESSION.get(full_url, timeout=timeout)
                resp.raise_for_status()
                data = resp.json()
                break
            except (requests.exceptions.JSONDecodeError, requests.exceptions.Timeout) as e:
                if attempt == retries - 1:
                    raise
                wait = 2 ** attempt
                logger.warning(f"Retry {attempt+1}/{retries}: {e}. Waiting {wait}s")
                time.sleep(wait)
            except Exception as e:
                if attempt == retries - 1:
                    raise
                wait = 2 ** attempt
                logger.warning(f"Retry {attempt+1}/{retries}: {e}. Waiting {wait}s")
                time.sleep(wait)

        if not data:
            break
        all_records.extend(data)
        batch_count += 1
        logger.debug(f"  offset={offset} fetched={len(data)} total={len(all_records)}")

        if batch_count >= save_every_batches:
            with open(partial_path, "w") as f:
                json.dump(all_records, f)
            logger.info(f"  Checkpoint: saved {len(all_records)} AMR records to {partial_path.name}")
            batch_count = 0

        content_range = resp.headers.get("Content-Range", "")
        if content_range:
            try:
                total = int(content_range.split("/")[-1])
                if offset + limit >= total:
                    break
            except ValueError:
                pass

        if len(data) < limit:
            break
        offset += limit
        time.sleep(0.2)

    return all_records


# ── Step 1: Genome metadata ───────────────────────────────────────────────────

def fetch_genome_metadata(cfg: dict) -> pd.DataFrame:
    logger.info("Fetching E. coli genome metadata from BV-BRC...")
    taxon = cfg["bvbrc"]["taxon_id"]
    qc = cfg["qc"]

    # Use a minimal API query (BV-BRC often 400s on complex RQL). Filter in Python.
    conditions = [
        f"eq(taxon_id,{taxon})",
        "eq(genome_quality,Good)",
    ]
    fields = ["genome_id", "genome_name", "genome_length", "contigs",
              "n50", "l50", "cds", "genome_status", "genome_quality", "assembly_accession"]

    records = bvbrc_query("genome", conditions, fields, limit=cfg["bvbrc"]["batch_size"])
    df = pd.DataFrame(records)
    if df.empty:
        logger.warning("No genomes returned from API")
        return df

    # Apply remaining QC filters in Python (only use columns that exist)
    for col in ["genome_length", "n50", "N50", "contigs"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    n50_col = "n50" if "n50" in df.columns else "N50" if "N50" in df.columns else None

    mask = pd.Series(True, index=df.index)
    if "genome_status" in df.columns:
        mask &= df["genome_status"].isin(["WGS", "Complete"])
    if "genome_length" in df.columns:
        mask &= (df["genome_length"] > qc["min_length"]) & (df["genome_length"] < qc["max_length"])
    if n50_col:
        mask &= df[n50_col] > qc["min_n50"]
    if "contigs" in df.columns:
        mask &= df["contigs"] < qc["max_contigs"]

    df = df[mask]
    logger.info(f"Found {len(df)} candidate genomes after QC filters")
    return df


# ── Step 2: AMR phenotypes ────────────────────────────────────────────────────

def fetch_amr_phenotypes(cfg: dict, out_amr_dir: Path) -> pd.DataFrame:
    logger.info("Fetching AMR phenotypes (resumable; progress saved every 100k records)...")
    taxon = cfg["bvbrc"]["taxon_id"]
    partial_path = out_amr_dir / "amr_phenotypes_partial.json"

    conditions = [
        f"eq(taxon_id,{taxon})",
        "in(resistant_phenotype,(Resistant,Susceptible))",
    ]
    fields = ["genome_id", "antibiotic", "resistant_phenotype",
              "laboratory_typing_method", "measurement_value", "measurement_units", "evidence"]

    records = bvbrc_query_resumable(
        "genome_amr", conditions, fields,
        limit=cfg["bvbrc"]["batch_size"],
        partial_path=partial_path,
        timeout=120,
        save_every_batches=20,
    )
    # Remove partial file on success so next run starts fresh (or uses full cache)
    if partial_path.exists():
        partial_path.unlink()
    df = pd.DataFrame(records)
    logger.info(f"Found {len(df)} AMR records covering {df['genome_id'].nunique()} genomes")
    return df


# ── Step 3: Select cohort ─────────────────────────────────────────────────────

def select_cohort(genome_df: pd.DataFrame, amr_df: pd.DataFrame, cfg: dict) -> list:
    target = cfg["cohort"]["target_genomes"]
    ab_interest = [a.lower().strip() for a in cfg["phenotypes"]["antibiotics_of_interest"]]

    amr = amr_df.copy()
    amr["antibiotic"] = amr["antibiotic"].str.lower().str.strip()

    if ab_interest:
        amr = amr[amr["antibiotic"].isin(ab_interest)]

    genomes_with_amr = set(amr["genome_id"].unique())
    available = set(genome_df["genome_id"].unique())
    intersection = genomes_with_amr & available
    logger.info(f"Genomes with both metadata and AMR data: {len(intersection)}")

    # Rank by number of antibiotics covered — maximises GWAS power
    counts = (
        amr[amr["genome_id"].isin(intersection)]
        .groupby("genome_id")["antibiotic"]
        .nunique()
        .sort_values(ascending=False)
    )
    selected = list(counts.head(target).index)
    logger.info(f"Selected {len(selected)} genomes for cohort")
    return selected


# ── Step 4: Download FASTAs ───────────────────────────────────────────────────

def download_fasta(genome_id: str, out_dir: Path, retries: int = 5) -> str:
    out_path = out_dir / f"{genome_id}.fna"
    if out_path.exists() and out_path.stat().st_size > 100_000:
        return "skipped"

    # BV-BRC FASTA endpoint — must set Accept header, not a JSON query
    url = f"https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,{genome_id})&select(sequence)&limit(5000,0)"
    headers = {"Accept": "application/dna+fasta"}

    for attempt in range(retries):
        try:
            resp = requests.get(url, headers=headers, timeout=120, stream=True)
            resp.raise_for_status()
            with open(out_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=65536):
                    f.write(chunk)
            size = out_path.stat().st_size
            if size < 1000:
                out_path.unlink(missing_ok=True)
                return "empty"
            return "downloaded"
        except Exception as e:
            if attempt == retries - 1:
                logger.error(f"Failed {genome_id}: {e}")
                out_path.unlink(missing_ok=True)
                return "failed"
            time.sleep(2 ** attempt)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    cfg = load_config()
    setup_logging(cfg)
    logger.info("=== Step 1: Download Data ===")

    out_genome_dir = ROOT / cfg["paths"]["raw_genomes"]
    out_amr_dir = ROOT / cfg["paths"]["raw_amr"]
    out_genome_dir.mkdir(parents=True, exist_ok=True)
    out_amr_dir.mkdir(parents=True, exist_ok=True)

    # Genome metadata
    genome_df = fetch_genome_metadata(cfg)
    genome_df.to_csv(out_amr_dir / "genome_metadata.tsv", sep="\t", index=False)

    # AMR phenotypes (cache to disk)
    amr_path = out_amr_dir / "amr_phenotypes.json"
    if amr_path.exists():
        logger.info(f"Loading cached AMR data from {amr_path}")
        with open(amr_path) as f:
            amr_df = pd.DataFrame(json.load(f))
    else:
        amr_df = fetch_amr_phenotypes(cfg, out_amr_dir)
        with open(amr_path, "w") as f:
            json.dump(amr_df.to_dict(orient="records"), f)

    # Cohort selection
    selected_ids = select_cohort(genome_df, amr_df, cfg)

    # Download FASTAs
    logger.info(f"Downloading {len(selected_ids)} genome FASTAs...")
    results = {}
    successful_ids = []

    for gid in tqdm(selected_ids, desc="Downloading genomes"):
        status = download_fasta(gid, out_genome_dir)
        results[status] = results.get(status, 0) + 1
        if status in ("downloaded", "skipped"):
            successful_ids.append(gid)
        time.sleep(0.2)

    logger.info(f"Download summary: {results}")

    # Save manifest
    manifest = genome_df[genome_df["genome_id"].isin(successful_ids)].copy()
    manifest_path = ROOT / "data" / "processed" / "genome_manifest.tsv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, sep="\t", index=False)

    print(f"\n=== Download Complete ===")
    print(f"Downloaded : {results.get('downloaded', 0)}")
    print(f"Skipped    : {results.get('skipped', 0)}")
    print(f"Failed     : {results.get('failed', 0)}")
    print(f"Empty      : {results.get('empty', 0)}")
    print(f"Manifest   : {manifest_path}")
    logger.info("=== Step 1 Complete ===")


if __name__ == "__main__":
    main()
