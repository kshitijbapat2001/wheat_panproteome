# Clean run as a SLURM job — `/scratch/pawsey1276/kbapat/Pan_prot_final/FINAL`

Two commands. Everything after setup runs inside one job.

```bash
cd /scratch/pawsey1276/kbapat/Pan_prot_final/FINAL
tar xzf PIPELINE_FINAL_JOB.tar.gz --strip-components=1

bash run_all.sh setup     # LOGIN NODE, ~20-40 min, needs internet
bash run_all.sh submit    # sbatch; the job does everything else
bash run_all.sh status    # queue + live log
bash run_all.sh finish    # post-processing on an existing build, if a job stopped early
bash run_all.sh collect   # bundle the outputs to send for the manuscript
```

## Your current build — finish it, do not re-run

The heavy stage completed successfully on 11 Aug:

```
outputs/pan_proteome_build_20260811_115422
  RCA   :  46 proteotypes (24 novel)   90 loci, 9 per cultivar, INVARIANT
  CPN60 : 141 proteotypes (89 novel)  345 loci, 32-36 per cultivar
  group-2 RCA clade: 10 proteotypes, monophyletic, support 0.99
```

The job then said `CRITICAL: no result directory produced` on a run that had in fact
worked. `config.yaml` carried **absolute** paths to the parent project, so the build was
written to `Pan_prot_final/outputs/` while the job looked in `FINAL/outputs/`. The paths
are now relative, and the job resolves the results directory from the config instead of
assuming. Nothing about the science was affected.

Recover that build without recomputing anything:

```bash
cd /scratch/pawsey1276/kbapat/Pan_prot_final/FINAL
bash finish_build.sh /scratch/pawsey1276/kbapat/Pan_prot_final/outputs/pan_proteome_build_20260811_115422
```

That runs audit, domain profiling, transcript cross-reference, audit again, figures,
tables, checks and the verification gate on the existing build — a few minutes, no
`sbatch`. `bash run_all.sh collect` then bundles it.

Note that because the old `config.yaml` pointed at the parent, this build sits in the
**old** outputs tree alongside earlier builds. That is exactly the situation that
produced the mixed-build bundle before, so let the gate run and check
`VERIFICATION.txt` before sending anything.

## What went wrong before that

`setup_env.sh` still demanded `build_manuscript.py`, which this build deliberately does
not ship, so setup aborted before creating the environment. Every later step then died
**silently**: `run_all.sh` used `set -e`, and its conda-activation helper returned
non-zero when conda did not exist, which killed the script before it printed anything.
Both are fixed — `run_all.sh` no longer uses `set -e` and every failure now names itself
and tells you which step to run.

## Step 1 — `bash run_all.sh setup` (login node)

Idempotent, safe to re-run. It:

1. installs **Miniforge** into `./miniforge3` (Setonix has no conda module, and HOME is
   too small for conda packages — all caches are redirected into the project directory);
2. creates `./envs/wheat_panproteome` from `environment.yml` — Python 3.10, DIAMOND,
   MAFFT, FastTree, HMMER, BioPython, pandas/numpy/scipy, matplotlib, pillow, openpyxl,
   python-docx;
3. downloads the ten Ensembl proteomes and the curated UniProt reference queries;
4. warms the UniProt/InterPro cache — **this is what outgroup rooting depends on**;
5. runs `prefetch_caches.py` for the Pfam HMMs, Swiss-Prot class seeds and the White
   et al. (2025) Zenodo tables.

Step 5 is new and it is what lets the job run offline. Compute nodes have no internet;
`domain_profiling.py` and `transcript_crossref.py` both download, so without warmed
caches they fail inside the job. `transcript_crossref.py` has also been patched to skip
the Zenodo API entirely when every file is already cached.

## Step 2 — `bash run_all.sh submit`

`submit_panproteome.sh` requests 64 CPUs, 115 GB, 24 h on `--partition=work`,
`--account=pawsey1276`, and runs, in order:

| | Step | Needs internet |
|---|---|---|
| 1 | `hpc_pipeline.py` — homology search, proteotypes, loci, alignment, trees, genome-wide clustering | no (cached) |
| 2 | `panproteome_audit.py` | no |
| 3 | `domain_profiling.py` — Pfam HMMs, class profiles, reciprocal best hits, CbbX control | no (prefetched) |
| 4 | `transcript_crossref.py` — White et al. cross-reference | no (prefetched) |
| 5 | `panproteome_audit.py` again, now that domain classes exist | no |
| 6 | `make_all_figures.py` — composite labels, 600 dpi | no |
| 7 | `build_tables.py` — Table 1–6, supplementary tables, `RESULTS_NUMBERS.json` | no |
| 8 | `presubmission_checks.py` | no |
| 9 | `verify_outputs.py` — **the gate** | no |

The job exits non-zero if verification fails, so a bad run cannot look like a good one
in the queue log.

**No manuscript prose is generated.** `build_manuscript.py` and `manuscript_text.py` are
not in this bundle; `build_tables.py` produces the data artefacts only.

## Step 3 — `bash run_all.sh collect`

Bundles `final/`, `audit/`, `domains/`, `transcripts/`, `checks/`, `tree_rooting.csv`,
the four core CSVs and `VERIFICATION.txt` into one tarball, and prints the `scp` line.

## The gate

`verify_outputs.py` fails on any defect that reached a draft during development:

- two copies of `proteotypes.csv` with different sizes (mixed builds)
- proteotype counts disagreeing across naming, domains and recovery tables
- outgroup rooting requested but midpoint achieved
- `Table7.csv` shorter than `candidate_table.csv`
- any blank cell in Table 1
- any figure below 300 dpi at 6.5-inch print width
- malformed composite identifiers

## Fixes carried in this build

1. **Figure DPI** — `hpc_pipeline.py` set no matplotlib DPI; its genome-wide figures
   saved at 100. Now 600 globally, explicit on every `savefig`.
2. **Panel letters** — the pipeline's pan-genome figure had no (a)/(b) although the
   caption cited them.
3. **Deterministic composite identifiers** — the representative cultivar was the first
   row among members tied at the best e-value, so reordering the input changed 13 of 175
   names. Ties now break on `(evalue, cultivar, source_id)`, in **both** copies of the
   function.
4. **Silent midpoint rooting** — a failed outgroup lookup fell back to midpoint with no
   error, contradicting the text. Now logs a loud failure and writes `tree_rooting.csv`.
5. **Blank Table 1 cells** — `reference_query_sets.csv` lives outside the build
   directory and the lookup missed it. Recursive fallback added, plus a warning.
6. **Figure numbering** — figures are copied into citation order, so
   `Figure5_pangenome_overview.png` becomes manuscript Figure 7 automatically.
7. **Offline transcript step** — skips the Zenodo API when the cache is complete.

## If something fails

```bash
bash run_all.sh status                     # queue + last 30 log lines
tail -100 logs/slurm_panproteome_<JID>.err
bash run_all.sh verify                     # re-run the gate alone
```

`setup` is safe to re-run; it reuses the environment and any files already downloaded.

## Still yours to decide

**Phenotype labels** are empty under `cultivars:` in `config.yaml`, so the association
framework stays reported as a future application. Set two or more cultivars to
`resistant` and two to `sensitive` and re-run to activate it.
