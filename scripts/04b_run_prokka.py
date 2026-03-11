#!/usr/bin/env python3
"""
04b_run_prokka.py  (SLOW PATH - optional)
Annotate each genome with Prokka in parallel.
Only needed if use_bvbrc_pa_matrix: false in config.yaml.
Requires: prokka installed (conda install -c bioconda prokka)
"""

import subprocess
from pathlib import Path

import yaml
from joblib import Parallel, delayed
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
    logger.add(log_dir / "04b_prokka.log", rotation="50 MB", level="DEBUG")


def load_passing_genomes(cfg: dict) -> list:
    path = ROOT / cfg["paths"]["passing_genomes"]
    return path.read_text().strip().splitlines()


def run_prokka(genome_id: str, fasta_path: Path, out_dir: Path, cpus: int, log_dir: Path) -> str:
    """Run Prokka for one genome. Returns 'done', 'skipped', or 'failed'."""
    genome_out = out_dir / genome_id
    gff_file = genome_out / f"{genome_id}.gff"

    if gff_file.exists() and gff_file.stat().st_size > 1000:
        return "skipped"

    genome_out.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"prokka_{genome_id}.log"

    cmd = [
        "prokka",
        "--outdir", str(genome_out),
        "--prefix", genome_id,
        "--genus", "Escherichia",
        "--species", "coli",
        "--kingdom", "Bacteria",
        "--cpus", str(cpus),
        "--quiet",
        "--force",
        str(fasta_path),
    ]

    with open(log_file, "w") as lf:
        result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        logger.error(f"Prokka failed for {genome_id} (exit {result.returncode})")
        return "failed"
    if not gff_file.exists():
        logger.error(f"Prokka output GFF missing for {genome_id}")
        return "failed"
    return "done"


def main():
    cfg = load_config()
    setup_logging(cfg)
    logger.info("=== Step 4b: Run Prokka (Slow Path) ===")

    if cfg["cohort"]["use_bvbrc_pa_matrix"]:
        logger.warning(
            "use_bvbrc_pa_matrix is True in config — Prokka not needed. "
            "Set it to False to use this script. Exiting."
        )
        return

    genome_ids = load_passing_genomes(cfg)
    genome_dir = ROOT / cfg["paths"]["raw_genomes"]
    prokka_dir = ROOT / cfg["paths"]["prokka_dir"]
    log_dir = ROOT / cfg["paths"]["logs"]
    prokka_dir.mkdir(parents=True, exist_ok=True)

    n_jobs = cfg["compute"]["prokka_jobs"]
    cpus_per_job = cfg["compute"]["prokka_cpus_per_job"]
    logger.info(f"Running {n_jobs} parallel Prokka jobs, {cpus_per_job} CPUs each")

    fasta_paths = [genome_dir / f"{gid}.fna" for gid in genome_ids]
    missing = [p for p in fasta_paths if not p.exists()]
    if missing:
        logger.warning(f"{len(missing)} FASTA files missing — skipping those genomes")

    valid = [(gid, p) for gid, p in zip(genome_ids, fasta_paths) if p.exists()]
    logger.info(f"Annotating {len(valid)} genomes...")

    results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(run_prokka)(gid, fasta, prokka_dir, cpus_per_job, log_dir)
        for gid, fasta in tqdm(valid, desc="Prokka")
    )

    from collections import Counter
    summary = Counter(results)
    logger.info(f"Prokka summary: {dict(summary)}")
    print(f"\nProkka done: {dict(summary)}")
    logger.info("=== Step 4b Complete ===")


if __name__ == "__main__":
    main()
