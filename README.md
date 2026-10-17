# wheat_panproteome

Pipeline that builds a cultivar-resolved pan-proteome for selected protein families from the annotated proteomes of bread wheat cultivars in Ensembl Plants, and writes the resulting data files, tables and figures. The configuration provided covers Rubisco activase (RCA) and chaperonin-60 (CPN60) in ten cultivars from Ensembl Plants release 63.

## Requirements

- Linux or macOS with `bash`
- conda (Miniforge is installed automatically if conda is not found)
- Internet access during setup
- Software environment defined in `environment.yml`: Python 3.10, DIAMOND, MAFFT, FastTree, HMMER, Biopython, pandas, NumPy, SciPy, Matplotlib, seaborn, PyYAML and openpyxl

## Usage

```bash
bash setup_env.sh               # environment, proteomes, reference sequences and caches (internet)
sbatch submit_panproteome.sh    # complete run on a SLURM cluster
bash submit_panproteome.sh      # or the same run on a workstation
```

`submit_panproteome.sh` runs offline once `setup_env.sh` has completed. Its `#SBATCH` header requests 64 CPUs, 115 GB of memory and 24 h; edit the account and partition lines for your cluster. To use an existing conda environment, set `CONDA_ENV_PREFIX` before running either script.

## Scripts

| Script | Function |
|---|---|
| `setup_env.sh` | Creates the conda environment; downloads the proteomes, UniProtKB reference sequences and annotation caches |
| `prefetch_caches.py` | Downloads the Pfam HMMs, Swiss-Prot seed sequences and independent annotation tables used offline |
| `submit_panproteome.sh` | Runs all analysis steps in order |
| `hpc_pipeline.py` | Homology search (DIAMOND), family assignment, proteoform collapsing, presence–absence matrix, gene loci, alignments (MAFFT), phylogenies (FastTree), physicochemical descriptors, per-column polymorphism, genome-wide clustering, annotated FASTA catalogues and SQLite database |
| `tree_tools.py` | Newick parsing, support summaries, collapsing of weakly supported nodes and rooting |
| `external_data.py` | UniProtKB and InterPro annotation; outgroup sequence retrieval |
| `panproteome_audit.py` | Gene-locus resolution, recruitment specificity, Cpn60 class assignment, divergence from the reference, prevalence at three clustering granularities and candidate list |
| `domain_profiling.py` | Pfam and class-specific HMM profiling with reciprocal best hits |
| `transcript_crossref.py` | Cross-reference of gene loci against an independent de novo annotation (Zenodo, doi:10.5281/zenodo.16964999) |
| `rca_isoforms.py` | Assignment of RCA proteoforms to the established RCA isoforms |
| `alignment_views.py` | Residue-level alignment views (PDF) |
| `make_nature_figures.py` | Figures (vector PDF and 600-dpi PNG) and source-data workbooks |
| `build_tables.py` | Tables, supplementary data files and the output bundle with an MD5 manifest |

## Configuration

All settings are in `config.yaml`: Ensembl release (`ensembl.release`), cultivar panel, UniProtKB queries for the reference sets, homology thresholds, whole-proteome clustering, alignment and rooting, outgroups, RCA isoform anchor loci and optional phenotype groups per cultivar (`group`). `make_nature_figures.py` draws the ten-cultivar panel defined in `config.yaml`.

## Outputs

Each run writes `outputs/pan_proteome_build_<timestamp>/`, linked as `outputs/latest`. The output bundle is written to `final/`:

```
final/
├── catalogue/              PanProteome_catalog.annotated.faa, WholeProteome_catalog.annotated.faa, wheat_pan_proteome.db
├── alignments_and_trees/   RCA and CPN60 alignments (.aln.fasta) and trees (.nwk, .supported.nwk)
├── reference/              reference_query_sets.csv
├── tables/                 Supplementary_Table_1-7.csv and Supplementary_Data_*.csv
├── supplementary_data/     Supplementary_Data_1-18 (.xlsx, .pdf, .zip)
├── figures/main/           Fig1-Fig4 (.pdf, .png)
├── figures/extended_data/  ExtendedDataFig1-3 (.pdf, .png)
├── source_data/            Source_Data_*.xlsx
└── MANIFEST.tsv            path, size and MD5 checksum of every file
```

Intermediate results remain in the build directory (`audit/`, `domains/`, `transcripts/`, `figures/`).

## Individual steps

Each step can be run on an existing build directory `BUILD`:

```bash
python panproteome_audit.py --inputs BUILD --outdir BUILD/audit
python domain_profiling.py --inputs BUILD --outdir BUILD/domains --cache inputs/hmm_cache
python transcript_crossref.py --inputs BUILD --outdir BUILD/transcripts --cache inputs/zenodo_cache
python rca_isoforms.py --inputs BUILD --outdir BUILD/audit --config config.yaml
python alignment_views.py --inputs BUILD --outdir BUILD/figures
python make_nature_figures.py --inputs BUILD --outdir BUILD/figures/nature --source-data BUILD/source_data
python build_tables.py --inputs BUILD --outdir BUILD/final --config config.yaml
```

`python hpc_pipeline.py --selftest` runs the offline unit tests.

## License

MIT; see `LICENSE`.

## Citation

See `CITATION.cff`.
