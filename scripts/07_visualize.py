#!/usr/bin/env python3
"""
07_visualize.py
Generate Manhattan plots, QQ plots, and resistance heatmap.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from loguru import logger
from scipy import stats

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.yaml"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def setup_logging(cfg):
    log_dir = ROOT / cfg["paths"]["logs"]
    log_dir.mkdir(parents=True, exist_ok=True)
    logger.add(log_dir / "07_visualize.log", rotation="10 MB", level="DEBUG")


def load_thresholds(cfg: dict) -> dict:
    thresh_path = ROOT / cfg["paths"]["gwas_inputs"] / "thresholds.txt"
    thresholds = {}
    if thresh_path.exists():
        for line in thresh_path.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                try:
                    thresholds[parts[0]] = float(parts[1])
                except ValueError:
                    thresholds[parts[0]] = parts[1]
    return thresholds


def load_gwas_results(antibiotic: str, gwas_results_dir: Path) -> pd.DataFrame | None:
    ab_safe = antibiotic.replace("/", "_").replace(" ", "_")
    results_path = gwas_results_dir / ab_safe / "gwas_results.tsv"
    if not results_path.exists():
        return None
    try:
        df = pd.read_csv(results_path, sep="\t", comment="#")
        pval_col = "lrt-pvalue" if "lrt-pvalue" in df.columns else "filter-pvalue"
        df["pvalue"] = df[pval_col]
        df["neg_log10_p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
        df["variant_index"] = range(len(df))
        return df
    except Exception as e:
        logger.error(f"Could not load results for {antibiotic}: {e}")
        return None


# ── Manhattan plot ─────────────────────────────────────────────────────────────

def plot_manhattan(df: pd.DataFrame, antibiotic: str, thresholds: dict, out_path: Path):
    fig, ax = plt.subplots(figsize=(14, 5))

    bonf = thresholds.get("bonferroni", 1e-5)
    sugg = thresholds.get("suggestive", 1e-5)
    bonf_log = -np.log10(bonf) if bonf > 0 else 8.0
    sugg_log = -np.log10(sugg) if sugg > 0 else 5.0

    # Color points by significance
    colors = []
    for _, row in df.iterrows():
        if row["neg_log10_p"] >= bonf_log:
            colors.append("#e74c3c")   # red: genome-wide significant
        elif row["neg_log10_p"] >= sugg_log:
            colors.append("#e67e22")   # orange: suggestive
        else:
            # Alternate grey shades for readability
            colors.append("#2c3e50" if row["variant_index"] % 2 == 0 else "#7f8c8d")

    ax.scatter(df["variant_index"], df["neg_log10_p"], c=colors, s=6, alpha=0.7, linewidths=0)

    # Threshold lines
    ax.axhline(y=bonf_log, color="#e74c3c", linestyle="--", linewidth=1.2, label=f"Bonferroni ({bonf:.1e})")
    ax.axhline(y=sugg_log, color="#e67e22", linestyle=":", linewidth=1.0, label=f"Suggestive ({sugg:.1e})")

    # Label top hits
    top_hits = df[df["neg_log10_p"] >= bonf_log].nlargest(10, "neg_log10_p")
    if not top_hits.empty:
        texts = []
        variant_col = "variant" if "variant" in df.columns else df.columns[0]
        for _, row in top_hits.iterrows():
            label = str(row[variant_col])[:20]
            t = ax.text(row["variant_index"], row["neg_log10_p"] + 0.1, label, fontsize=6)
            texts.append(t)
        if HAS_ADJUSTTEXT and texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))

    ax.set_xlabel("Gene cluster (index)", fontsize=11)
    ax.set_ylabel(r"$-\log_{10}(p)$", fontsize=11)
    ax.set_title(f"GWAS Manhattan Plot — {antibiotic.title()}", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")

    # Legend for colors
    patches = [
        mpatches.Patch(color="#e74c3c", label="Significant"),
        mpatches.Patch(color="#e67e22", label="Suggestive"),
        mpatches.Patch(color="#7f8c8d", label="Not significant"),
    ]
    ax.legend(handles=patches, fontsize=9, loc="upper right")

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Manhattan plot saved: {out_path}")


# ── QQ plot ───────────────────────────────────────────────────────────────────

def plot_qq(df: pd.DataFrame, antibiotic: str, out_path: Path):
    pvals = df["pvalue"].dropna()
    pvals = pvals[pvals > 0].sort_values()
    n = len(pvals)

    observed = -np.log10(pvals.values)
    expected = -np.log10(np.linspace(1 / n, 1, n))

    # Genomic inflation factor λ
    chi2_obs = stats.chi2.ppf(1 - pvals.values, df=1)
    lambda_gc = np.median(chi2_obs) / stats.chi2.ppf(0.5, df=1)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=8, alpha=0.6, color="#2980b9", linewidths=0)
    max_val = max(observed.max(), expected.max()) * 1.05
    ax.plot([0, max_val], [0, max_val], color="#e74c3c", linewidth=1.5, linestyle="--", label="Expected")
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$", fontsize=11)
    ax.set_ylabel(r"Observed $-\log_{10}(p)$", fontsize=11)
    ax.set_title(f"QQ Plot — {antibiotic.title()}\n$\\lambda_{{GC}}$ = {lambda_gc:.3f}", fontsize=12)
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"QQ plot saved: {out_path}")


# ── Resistance heatmap ────────────────────────────────────────────────────────

def plot_resistance_heatmap(pheno_matrix: pd.DataFrame, out_path: Path):
    if pheno_matrix.empty:
        logger.warning("Phenotype matrix is empty — skipping heatmap")
        return

    # Drop genomes with all-NaN phenotypes
    df = pheno_matrix.dropna(how="all").dropna(axis=1, how="all")
    if df.empty:
        return

    # Fill NaN with -1 for visualization
    df_plot = df.fillna(-1)

    cmap = plt.cm.colors.ListedColormap(["#bdc3c7", "#27ae60", "#e74c3c"])
    bounds = [-1.5, -0.5, 0.5, 1.5]
    norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)

    fig_height = max(8, len(df_plot) * 0.08)
    fig_width = max(8, len(df_plot.columns) * 1.2)
    fig, ax = plt.subplots(figsize=(min(fig_width, 20), min(fig_height, 20)))

    im = ax.imshow(df_plot.values, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(df_plot.columns)))
    ax.set_xticklabels(df_plot.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks([])
    ax.set_xlabel("Antibiotic", fontsize=12)
    ax.set_ylabel(f"Genomes (n={len(df_plot)})", fontsize=12)
    ax.set_title("E. coli Antibiotic Resistance Phenotypes", fontsize=13, fontweight="bold")

    # Legend
    legend_elements = [
        mpatches.Patch(color="#e74c3c", label="Resistant"),
        mpatches.Patch(color="#27ae60", label="Susceptible"),
        mpatches.Patch(color="#bdc3c7", label="Missing"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1.15, 1), fontsize=10)

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Resistance heatmap saved: {out_path}")


# ── Summary table ─────────────────────────────────────────────────────────────

def generate_summary_table(gwas_results_dir: Path, antibiotics: list, thresholds: dict) -> pd.DataFrame:
    bonf = thresholds.get("bonferroni", 1e-5)
    rows = []
    for ab in antibiotics:
        ab_safe = ab.replace("/", "_").replace(" ", "_")
        sig_path = gwas_results_dir / ab_safe / "significant_hits.tsv"
        if sig_path.exists():
            sig_df = pd.read_csv(sig_path, sep="\t")
            for _, row in sig_df.iterrows():
                pval_col = "lrt-pvalue" if "lrt-pvalue" in sig_df.columns else "filter-pvalue"
                rows.append({
                    "antibiotic": ab,
                    "variant": row.get("variant", ""),
                    "p_value": row.get(pval_col, np.nan),
                    "beta": row.get("beta", np.nan),
                    "allele_freq": row.get("af", np.nan),
                })
    return pd.DataFrame(rows)


def main():
    cfg = load_config()
    setup_logging(cfg)
    logger.info("=== Step 7: Visualize Results ===")

    figures_dir = ROOT / cfg["paths"]["figures"]
    figures_dir.mkdir(parents=True, exist_ok=True)
    gwas_results_dir = ROOT / cfg["paths"]["gwas_results"]
    pheno_path = ROOT / cfg["paths"]["phenotype_matrix"]
    thresholds = load_thresholds(cfg)

    if not pheno_path.exists():
        logger.error(f"Phenotype matrix not found: {pheno_path}. Run 03_prepare_phenotypes.py first.")
        return

    pheno_df = pd.read_csv(pheno_path, sep="\t", index_col=0)
    antibiotics = list(pheno_df.columns)

    # 1. Resistance heatmap (always possible after step 3)
    plot_resistance_heatmap(pheno_df, figures_dir / "heatmap_resistance.png")

    # 2. Per-antibiotic Manhattan + QQ plots
    any_gwas = False
    for ab in antibiotics:
        df = load_gwas_results(ab, gwas_results_dir)
        if df is None or df.empty:
            logger.warning(f"No GWAS results for {ab} — skipping plots")
            continue
        any_gwas = True
        ab_safe = ab.replace("/", "_").replace(" ", "_")
        plot_manhattan(df, ab, thresholds, figures_dir / f"manhattan_{ab_safe}.png")
        plot_qq(df, ab, figures_dir / f"qq_{ab_safe}.png")

    # 3. Summary table
    if any_gwas:
        summary = generate_summary_table(gwas_results_dir, antibiotics, thresholds)
        if not summary.empty:
            summary_path = ROOT / "results" / "significant_hits_summary.tsv"
            summary.to_csv(summary_path, sep="\t", index=False)
            logger.info(f"Summary table saved: {summary_path}")

            print("\n=== Top Significant Associations ===")
            print(summary.sort_values("p_value").head(20).to_string(index=False))
        else:
            logger.info("No significant hits to summarize.")
    else:
        logger.info("No GWAS results found — only heatmap was generated.")
        logger.info("Run 06_run_gwas.py first to generate Manhattan and QQ plots.")

    logger.info("=== Step 7 Complete ===")
    print(f"\nFigures saved to: {figures_dir}")


if __name__ == "__main__":
    main()
