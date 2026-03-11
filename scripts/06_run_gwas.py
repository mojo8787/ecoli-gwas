#!/usr/bin/env python3
"""
06_run_gwas.py
Run pyseer GWAS for each antibiotic with LMM population structure correction.

Steps:
  1. Build Mash distance matrix (population structure)
  2. Convert to similarity matrix for pyseer
  3. Count patterns (for Bonferroni correction)
  4. Run pyseer LMM per antibiotic
  5. Apply Bonferroni threshold and save significant hits
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
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
    logger.add(log_dir / "06_gwas.log", rotation="50 MB", level="DEBUG")


def load_passing_genomes(cfg: dict) -> list:
    path = ROOT / cfg["paths"]["passing_genomes"]
    return path.read_text().strip().splitlines()


def find_tool(name: str) -> str | None:
    """Find a tool, preferring the venv bin directory alongside this Python."""
    venv_bin = Path(sys.executable).parent / name
    if venv_bin.exists():
        return str(venv_bin)
    result = subprocess.run(["which", name], capture_output=True, text=True)
    if result.returncode == 0:
        return result.stdout.strip()
    return None


def check_tool(name: str) -> bool:
    return find_tool(name) is not None


# ── Step 1: Mash distances ────────────────────────────────────────────────────

def build_mash_distances(cfg: dict) -> Path:
    """Build a Mash distance matrix. Returns path to distances TSV."""
    genome_dir = ROOT / cfg["paths"]["raw_genomes"]
    struct_dir = ROOT / cfg["paths"]["gwas_inputs"] / "population_structure"
    struct_dir.mkdir(parents=True, exist_ok=True)

    sketch_path = struct_dir / "sketches"
    dist_path = struct_dir / "mash_distances.tsv"
    sim_path = struct_dir / "mash_similarity.tsv"

    if sim_path.exists():
        logger.info(f"Mash similarity matrix already exists: {sim_path}")
        return sim_path

    fasta_files = sorted(genome_dir.glob("*.fna"))
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files in {genome_dir}")

    # Sketch
    logger.info(f"Running mash sketch on {len(fasta_files)} genomes...")
    cmd_sketch = ["mash", "sketch", "-s", "10000", "-o", str(sketch_path)] + [str(f) for f in fasta_files]
    result = subprocess.run(cmd_sketch, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"mash sketch failed: {result.stderr}")
        raise RuntimeError("mash sketch failed")

    # Distance matrix
    logger.info("Computing Mash pairwise distances...")
    sketch_msh = Path(str(sketch_path) + ".msh")
    with open(dist_path, "w") as f:
        result = subprocess.run(
            ["mash", "dist", "-t", str(sketch_msh), str(sketch_msh)],
            stdout=f, stderr=subprocess.PIPE, text=True
        )
    if result.returncode != 0:
        logger.error(f"mash dist failed: {result.stderr}")
        raise RuntimeError("mash dist failed")

    # Convert distance → similarity and strip file paths from genome IDs
    logger.info("Converting Mash distances to similarity matrix...")
    dist_df = pd.read_csv(dist_path, sep="\t", index_col=0)
    # Strip path prefix and .fna suffix from genome IDs
    dist_df.index = [Path(idx).stem for idx in dist_df.index]
    dist_df.columns = [Path(c).stem for c in dist_df.columns]
    sim_df = 1.0 - dist_df
    sim_df.to_csv(sim_path, sep="\t")
    logger.info(f"Similarity matrix saved: {sim_path}")
    return sim_path


def build_python_similarity(cfg: dict) -> Path:
    """
    Fallback: build a pairwise similarity matrix using Biopython-based
    Jaccard distance on the gene PA matrix (used when mash is not installed).
    """
    from sklearn.metrics import pairwise_distances
    logger.info("mash not found — building similarity matrix from gene PA matrix (fallback)")

    pa_path = ROOT / cfg["paths"]["gene_pa_matrix"]
    pa = pd.read_csv(pa_path, sep="\t", index_col=0)

    # Jaccard similarity on binary presence/absence
    dist = pairwise_distances(pa.values, metric="jaccard", n_jobs=cfg["compute"]["python_jobs"])
    sim = 1.0 - dist
    sim_df = pd.DataFrame(sim, index=pa.index, columns=pa.index)

    struct_dir = ROOT / cfg["paths"]["gwas_inputs"] / "population_structure"
    struct_dir.mkdir(parents=True, exist_ok=True)
    sim_path = struct_dir / "mash_similarity.tsv"
    sim_df.to_csv(sim_path, sep="\t")
    logger.info(f"Fallback similarity matrix saved: {sim_path}")
    return sim_path


# ── Step 2: Pattern count for Bonferroni threshold ───────────────────────────

def count_patterns(pa_path: Path, pheno_path: Path, gwas_inputs: Path) -> int:
    """Run pyseer count_patterns to get effective number of tests."""
    pattern_file = gwas_inputs / "pattern_count.txt"

    if pattern_file.exists():
        with open(pattern_file) as f:
            for line in f:
                if "patterns" in line.lower():
                    try:
                        return int(line.split()[-1])
                    except ValueError:
                        pass

    logger.info("Counting unique patterns for Bonferroni correction...")
    cmd = [
        sys.executable, "-m", "pyseer.scripts.count_patterns",
        "--phenotypes", str(pheno_path),
        "--pres", str(pa_path),
    ]

    with open(pattern_file, "w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logger.warning(f"count_patterns failed: {result.stderr}. Using total feature count.")
        # pa_path is Rtab: rows=genes, columns=samples → shape[0] = number of genes
        pa = pd.read_csv(pa_path, sep="\t", index_col=0)
        return pa.shape[0]

    with open(pattern_file) as f:
        for line in f:
            if "patterns" in line.lower():
                try:
                    return int(line.split()[-1])
                except ValueError:
                    pass

    pa = pd.read_csv(pa_path, sep="\t", index_col=0)
    return pa.shape[0]


# ── Step 3: Run pyseer per antibiotic ────────────────────────────────────────

def run_pyseer_for_antibiotic(
    antibiotic: str,
    pheno_path: Path,
    pa_path: Path,
    sim_path: Path,
    out_dir: Path,
    cfg: dict,
) -> Path:
    """Run pyseer LMM for one antibiotic. Returns path to results TSV."""
    ab_dir = out_dir / antibiotic.replace("/", "_").replace(" ", "_")
    ab_dir.mkdir(parents=True, exist_ok=True)
    results_path = ab_dir / "gwas_results.tsv"
    patterns_path = ab_dir / "patterns.txt"

    if results_path.exists() and results_path.stat().st_size > 100:
        logger.info(f"  {antibiotic}: results already exist, skipping")
        return results_path

    logger.info(f"  Running pyseer for: {antibiotic}")

    use_lmm = cfg["gwas"]["lmm"]
    cpus = cfg["compute"]["pyseer_cpus"]

    pyseer_bin = find_tool("pyseer")
    cmd = [
        pyseer_bin,
        "--lmm" if use_lmm else "--linear",
        "--phenotypes", str(pheno_path),
        "--phenotype-column", antibiotic,
        "--pres", str(pa_path),
        "--similarity", str(sim_path),
        "--output-patterns", str(patterns_path),
        "--cpu", str(cpus),
        "--min-af", str(cfg["gwas"]["min_af"]),
        "--max-af", str(cfg["gwas"]["max_af"]),
    ]

    log_file = ROOT / cfg["paths"]["logs"] / f"pyseer_{antibiotic.replace('/', '_')}.log"
    with open(results_path, "w") as out_f, open(log_file, "w") as log_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=log_f)

    if result.returncode != 0:
        logger.error(f"  pyseer failed for {antibiotic}. Check {log_file}")
        return None

    logger.info(f"  {antibiotic}: GWAS complete → {results_path}")
    return results_path


def extract_significant_hits(results_path: Path, threshold: float, antibiotic: str) -> pd.DataFrame:
    """Load pyseer results and return rows below the significance threshold."""
    try:
        df = pd.read_csv(results_path, sep="\t", comment="#")
    except Exception as e:
        logger.error(f"Could not read {results_path}: {e}")
        return pd.DataFrame()

    # pyseer output columns: variant, af, filter-pvalue, lrt-pvalue, beta, ...
    pval_col = "lrt-pvalue" if "lrt-pvalue" in df.columns else "filter-pvalue"
    df["antibiotic"] = antibiotic
    df["neg_log10_p"] = -np.log10(df[pval_col].clip(lower=1e-300))

    significant = df[df[pval_col] < threshold].copy()
    significant = significant.sort_values(pval_col)
    return significant


def main():
    cfg = load_config()
    setup_logging(cfg)
    logger.info("=== Step 6: Run GWAS ===")

    pa_matrix_path = ROOT / cfg["paths"]["gene_pa_matrix"]
    pheno_path = ROOT / cfg["paths"]["phenotype_matrix"]
    gwas_results_dir = ROOT / cfg["paths"]["gwas_results"]
    gwas_results_dir.mkdir(parents=True, exist_ok=True)
    gwas_inputs = ROOT / cfg["paths"]["gwas_inputs"]

    if not pa_matrix_path.exists():
        logger.error(f"Gene PA matrix not found: {pa_matrix_path}. Run 04_build_pa_matrix.py first.")
        return

    # pyseer --pres requires Rtab format: rows=genes, columns=samples, header "Gene\tsample1\tsample2..."
    # gene_pa_matrix.tsv is transposed (rows=samples, columns=genes), so we transpose it.
    pa_path = gwas_inputs / "gene_pa_rtab.tsv"
    if not pa_path.exists():
        logger.info("Transposing gene PA matrix to Rtab format for pyseer...")
        pa_df = pd.read_csv(pa_matrix_path, sep="\t", index_col=0)
        rtab = pa_df.T
        rtab.index.name = "Gene"
        rtab.to_csv(pa_path, sep="\t")
        logger.info(f"Rtab file written: {pa_path}")
    if not pheno_path.exists():
        logger.error(f"Phenotype matrix not found: {pheno_path}. Run 03_prepare_phenotypes.py first.")
        return

    pheno_df = pd.read_csv(pheno_path, sep="\t", index_col=0)
    antibiotics = list(pheno_df.columns)
    logger.info(f"Antibiotics to test: {antibiotics}")

    # Population structure
    if check_tool("mash"):
        sim_path = build_mash_distances(cfg)
    else:
        logger.warning("mash not found — using fallback similarity matrix from gene PA matrix")
        sim_path = build_python_similarity(cfg)

    # Bonferroni threshold (count_patterns also needs the Rtab file)
    n_patterns = count_patterns(pa_path, pheno_path, gwas_inputs)
    alpha = cfg["gwas"]["significance_alpha"]
    threshold = alpha / max(n_patterns, 1)
    suggestive = cfg["gwas"]["suggestive_threshold"]
    logger.info(f"Patterns: {n_patterns}  |  Bonferroni threshold: {threshold:.2e}")

    # Save threshold for visualizer
    (gwas_inputs / "thresholds.txt").write_text(
        f"bonferroni\t{threshold}\nsuggestive\t{suggestive}\npatterns\t{n_patterns}\n"
    )

    # Per-antibiotic GWAS
    if not check_tool("pyseer"):
        logger.error("pyseer not found. Install with: pip install pyseer")
        return

    all_significant = []
    for ab in antibiotics:
        # Drop genomes with missing phenotype for this antibiotic
        col = pheno_df[ab].dropna()
        if len(col) < 40:
            logger.warning(f"  {ab}: only {len(col)} samples with phenotype — skipping")
            continue

        res_path = run_pyseer_for_antibiotic(
            antibiotic=ab,
            pheno_path=pheno_path,
            pa_path=pa_path,
            sim_path=sim_path,
            out_dir=gwas_results_dir,
            cfg=cfg,
        )
        if res_path:
            sig = extract_significant_hits(res_path, threshold, ab)
            if not sig.empty:
                ab_dir = gwas_results_dir / ab.replace("/", "_").replace(" ", "_")
                sig.to_csv(ab_dir / "significant_hits.tsv", sep="\t", index=False)
                all_significant.append(sig)
                logger.info(f"  {ab}: {len(sig)} significant hits")
            else:
                logger.info(f"  {ab}: no significant hits at Bonferroni threshold")

    # Combined significant hits
    if all_significant:
        combined = pd.concat(all_significant, ignore_index=True)
        combined_path = gwas_results_dir / "all_significant_hits.tsv"
        combined.to_csv(combined_path, sep="\t", index=False)
        logger.info(f"Combined significant hits saved: {combined_path}")
        print(f"\nTotal significant associations: {len(combined)}")
        print(combined[["antibiotic", "variant", "lrt-pvalue", "beta"]].head(20).to_string(index=False))
    else:
        logger.warning("No significant hits found across any antibiotic.")

    logger.info("=== Step 6 Complete ===")


if __name__ == "__main__":
    main()
