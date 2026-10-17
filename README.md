# Wheat RCA/Cpn60 pan-proteome — code and data release (v2.0)

Code and data for **"Pan-proteome catalogue of Rubisco activase and chaperonins in bread wheat"** (K. Bapat, A. P. Scafaro, A. Pratap and N. L. Taylor; submitted to *Nature Plants* as a Resource).

Dataset DOI (Zenodo): https://doi.org/10.5281/zenodo.21886022 (concept DOI for all versions: https://doi.org/10.5281/zenodo.21886021)
Code: https://github.com/kshitijbapat2001/wheat_panproteome

The release contains (i) the complete analysis pipeline and (ii) every derived dataset, figure, table, Supplementary Data file and Source Data file of the manuscript.

## Contents

| Path | Contents |
|---|---|
| `code/` | Pipeline (Python 3, Bash, SLURM): homology search, proteoform collapsing, gene-locus resolution, audit, domain profiling, independent-annotation cross-reference, RCA isoform assignment, alignments, phylogenies, figures and tables |
| `docs/` | `RUNBOOK.md` (running the pipeline) and `RERUN.md` (clean re-run from an empty environment) |
| `data/catalogue/` | `PanProteome_catalog.annotated.faa` (187 RCA and Cpn60 proteoforms, composite identifiers in the headers), `WholeProteome_catalog.annotated.faa` (131,282 genome-wide proteoform clusters), `wheat_pan_proteome.db` (SQLite database of the catalogue) |
| `data/alignments_and_trees/` | MAFFT alignments (FASTA) and FastTree trees (Newick): raw (`*.nwk`) and support-collapsed, outgroup-rooted (`*.supported.nwk`) |
| `data/reference/` | The 64 curated UniProtKB reference queries (28 RCA, 36 Cpn60) |
| `data/tables/` | Plain CSV versions of Supplementary Tables 1–7 and of every tabular Supplementary Data item |
| `data/supplementary_data/` | Supplementary Data 1–18 exactly as submitted (Excel workbooks with a `Legend` sheet describing every column; PDFs; ZIP) |
| `data/figures/` | Fig. 1–4 and Extended Data Figs 1-3 (vector PDF with embedded fonts and 600-dpi PNG) |
| `data/source_data/` | Source Data for every figure (Excel, one workbook per figure) |
| `MANIFEST.tsv` | Every file with size and MD5 checksum |


## Running the pipeline

```bash
conda env create -f code/environment.yml        # or: bash code/setup_env.sh
conda activate wheat_panproteome
cd code
bash setup_env.sh                                 # proteomes, reference queries, caches (needs internet)
sbatch submit_panproteome.sh                      # full build on a SLURM cluster
bash run_all.sh status
```

The build writes `outputs/pan_proteome_build_<timestamp>/`. `verify_outputs.py` is the consistency gate. The figures in this release are produced from a build directory with

```bash
python make_nature_figures.py --inputs <build> --outdir figures --source-data source_data
python alignment_views.py --inputs <build> --outdir supplementary   # Supplementary Data 9 and 12
```

Input proteomes: Ensembl Plants release 63 (https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants) for Chinese Spring, ArinaLrFor, Jagger, Julius, Lancer, Landmark, Mace, SY Mattis, Norin61 and CDC Stanley. Independent annotation: White et al. (2025), Zenodo https://doi.org/10.5281/zenodo.16964999.

## Licences

Code (`code/`): MIT licence (`LICENSE`). Data, figures and tables (`data/`): Creative Commons Attribution 4.0 International (`LICENSE-DATA`).

## Citation

Bapat, K., Scafaro, A. P., Pratap, A. & Taylor, N. L. Pan-proteome catalogue of Rubisco activase and chaperonins in bread wheat. *Nature Plants* (submitted). Data and code: Zenodo https://doi.org/10.5281/zenodo.21886022.

Contact: Agyeya Pratap (agyeya.pratap@uwa.edu.au) and Nicolas L. Taylor (nicolas.taylor@uwa.edu.au).
