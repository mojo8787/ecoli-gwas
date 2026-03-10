# E. coli GWAS Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

A fully automated pipeline for performing Genome-Wide Association Studies (GWAS) on *Escherichia coli* to identify genomic loci (genes) associated with antibiotic resistance.

## What it does

This pipeline downloads publicly available *E. coli* genomes and AMR (antimicrobial resistance) phenotypes from the [BV-BRC database](https://www.bv-brc.org/), runs quality control, builds a pangenome presence/absence matrix, and uses **pyseer** with a Linear Mixed Model (LMM) to correct for population structure while identifying resistance-associated genes. Results are visualised as Manhattan plots, QQ plots, and a resistance heatmap.

### Antibiotics tested (by default)
ciprofloxacin, ampicillin, tetracycline, trimethoprim, ceftriaxone, meropenem, gentamicin, chloramphenicol

---

## Pipeline steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_download_data.py` | Download up to 200 *E. coli* genomes + AMR phenotypes from BV-BRC API |
| 2 | `02_qc_genomes.py` | Filter genomes by length, contig count, N50, and GC content |
| 3 | `03_prepare_phenotypes.py` | Build a binary (R/S) phenotype matrix per antibiotic |
| 4 | `04_build_pa_matrix.py` | Build gene presence/absence matrix via BV-BRC plfam API (fast path) or Prokka + Panaroo (slow path) |
| 5 | `06_run_gwas.py` | Compute Mash distances for population structure, run pyseer LMM per antibiotic, apply Bonferroni correction |
| 6 | `07_visualize.py` | Generate Manhattan plots, QQ plots, and resistance heatmap |

---

## Requirements

### System tools (must be on PATH)
- `mash` — for computing genomic distances
- `prokka` — genome annotation (only needed for slow path)
- `panaroo` — pangenome analysis (only needed for slow path)

### Python packages
```
pip install -r requirements.txt
```

Or use the provided conda environment:
```
conda env create -f environment.yml
conda activate ecoli_gwas
```

---

## Quick start

```bash
# 1. Clone / navigate to the project
cd ecoli_gwas

# 2. Create and activate a virtual environment
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# 3. (Optional) Edit config/config.yaml to adjust cohort size, QC thresholds, antibiotics, etc.

# 4. Run the full pipeline
bash run_pipeline.sh
```

---

## Configuration

All parameters are in `config/config.yaml`:

| Section | Key | Default | Description |
|---------|-----|---------|-------------|
| `cohort` | `target_genomes` | 200 | Max genomes to download |
| `cohort` | `use_bvbrc_pa_matrix` | `true` | Fast path (BV-BRC) vs slow path (Prokka+Panaroo) |
| `qc` | `min_length` / `max_length` | 4–6.5 Mbp | Genome size filter |
| `qc` | `max_contigs` | 500 | Assembly fragmentation filter |
| `qc` | `min_n50` | 10 000 bp | Assembly quality filter |
| `gwas` | `lmm` | `true` | Use LMM for population structure correction |
| `gwas` | `significance_alpha` | 0.05 | Bonferroni significance level |
| `phenotypes` | `min_resistant` / `min_susceptible` | 20 | Minimum class size per antibiotic |

---

## Outputs

```
results/
  figures/
    manhattan_<antibiotic>.png   # Manhattan plot per antibiotic
    qq_<antibiotic>.png          # QQ plot per antibiotic
    resistance_heatmap.png       # Heatmap of resistance prevalence
  significant_hits_summary.tsv  # All significant gene associations

data/gwas/results/
  <antibiotic>_gwas.tsv         # Full pyseer output per antibiotic
  <antibiotic>_significant.tsv  # Significant hits per antibiotic

logs/                            # Per-step log files
```

---

## Project structure

```
ecoli_gwas/
├── config/
│   └── config.yaml          # All pipeline parameters
├── data/
│   ├── raw/                 # Downloaded genomes and AMR data
│   ├── processed/           # QC-filtered genomes, phenotypes, pangenome
│   └── gwas/                # GWAS inputs and results
├── scripts/
│   ├── 01_download_data.py
│   ├── 02_qc_genomes.py
│   ├── 03_prepare_phenotypes.py
│   ├── 04_build_pa_matrix.py
│   ├── 04b_run_prokka.py    # Optional: Prokka annotation (slow path)
│   ├── 05_run_panaroo.py    # Optional: Panaroo pangenome (slow path)
│   ├── 06_run_gwas.py
│   └── 07_visualize.py
├── notebooks/               # Exploratory Jupyter notebooks
├── results/                 # Final figures and summary tables
├── logs/                    # Pipeline logs
├── run_pipeline.sh          # One-command pipeline runner
├── requirements.txt
└── environment.yml
```

---

## Data source

Genomes and AMR phenotypes are fetched from the **BV-BRC** (Bacterial and Viral Bioinformatics Resource Center) public REST API using taxon ID `562` (*E. coli*).

---

## Citation

If you use this pipeline in your research, please cite the software and the Zenodo archive (DOI will appear after your first release):

- **Software**: See [CITATION.cff](CITATION.cff) or the "Cite this repository" widget on GitHub.
- **Zenodo**: After connecting this repo to [Zenodo](https://zenodo.org) and creating a release, a DOI will be assigned. Update the badge above with your DOI (replace `XXXXXXX`).

---

## Method notes

- **Population structure correction**: Mash sketch distances are converted to a similarity matrix and used as the kinship matrix in pyseer's LMM mode, reducing false positives caused by phylogenetic relatedness.
- **Multiple testing**: Bonferroni correction is applied using the number of unique k-mer/gene patterns rather than raw variant count, as recommended by pyseer.
- **Fast vs slow pangenome path**: The fast path queries BV-BRC's pre-computed protein family (plfam) presence/absence, avoiding local Prokka+Panaroo runs. Set `use_bvbrc_pa_matrix: false` in the config to use Prokka+Panaroo instead.
