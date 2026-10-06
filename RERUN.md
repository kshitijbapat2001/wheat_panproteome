# Complete re-run from a clean environment

**Everything below runs in `/scratch/pawsey1276/kbapat/Pan_prot_final/FINAL`.**

Every fix found during development is baked into this bundle. Tables and figures are
numbered in the order they are first cited in the manuscript (v2.1): main Tables 1-7,
Supplementary Tables S1-S16 and Supplementary Figures S1-S6. The order is defined once, in
`numbering.py`, and every output file name follows it. Supplementary tables are written in
publication format with the two required header lines and a per-column description block.

**Already have a verified build?** You do not need to rerun anything heavy. Run
`bash apply_revision_v21.sh` (about two minutes) and skip to section 5.

---

## 1. Put the bundle in place

```bash
cd /scratch/pawsey1276/kbapat/Pan_prot_final/FINAL
tar xzf PIPELINE_RERUN.tar.gz --strip-components=1
ls *.py *.sh *.yaml *.yml | wc -l          # 29
```

Overwrites the code only. `inputs/` and `outputs/` are untouched at this stage.

## 2. See what the reset will remove — change nothing yet

```bash
bash reset_and_run.sh --plan
```

Prints, with sizes, exactly what will be deleted and what will be kept. Read it before
continuing. Nothing is modified by `--plan`.

| Removed | Kept |
|---|---|
| `envs/` `miniforge3/` `outputs/` `logs/` `__pycache__/` | every `.py` `.sh` `.yaml` `.yml` `.md` in this directory |
| `github_stage/` `zenodo_upload/` `review_*/` | `inputs/` — proteomes and warmed caches |
| `*.tar.gz` `*.zip` `*.bak` `*.pre_v3` `*.orig.png` | `archive/` — see below |

`inputs/` is kept because it holds the ten downloaded proteomes and the UniProt, Pfam and
Zenodo caches. Deleting them costs an hour of downloads and changes nothing about the
result. Use `--purge-inputs` if you want a genuinely bare start.

**Any existing `outputs/` is archived to `archive/outputs_<stamp>.tar.gz` before removal.**
Nothing is destroyed; the previous build stays recoverable.

## 3. Run it

```bash
bash reset_and_run.sh
```

Type `RESET` at the prompt. The script then:

1. archives and removes the directories above
2. runs `setup_env.sh` — installs Miniforge into `./miniforge3`, builds
   `./envs/wheat_panproteome` from `environment.yml`, downloads the ten proteomes and the
   curated UniProt reference queries, warms the UniProt/InterPro cache, and pre-fetches
   the Pfam HMMs, Swiss-Prot seeds and the White et al. Zenodo tables
3. verifies `inputs/hmm_cache` and `inputs/zenodo_cache` exist, and **stops** if they do
   not, since compute nodes have no internet and the domain and transcript steps would
   fail inside the job
4. submits `submit_panproteome.sh` and waits
5. runs `verify_outputs.py`
6. writes the publication-formatted supplementary tables

Steps 1–2 need a login node. Expect 30–60 minutes for setup, then the job itself.

If you would rather drive it yourself, the same sequence is:

```bash
bash reset_and_run.sh --plan          # inspect
bash reset_and_run.sh --yes           # no prompt
# or, manually:
bash setup_env.sh
sbatch submit_panproteome.sh
bash run_all.sh status
```

## 4. What the job runs

| | Step | Internet |
|---|---|---|
| 1 | `hpc_pipeline.py` — homology search, proteotypes, gene loci, alignments, trees, genome-wide clustering | cached |
| 2 | `panproteome_audit.py` | no |
| 3 | `domain_profiling.py` — Pfam HMMs, class profiles, reciprocal best hits, CbbX control | prefetched |
| 4 | `transcript_crossref.py` — White et al. cross-reference | prefetched |
| 5 | `panproteome_audit.py` again, now that domain classes exist | no |
| 5b | `rca_isoforms.py` — RCA proteotypes to the nine established isoforms (Table 5, S7) | no |
| 5c | `alignment_views.py` — residue-level alignment figures (S3, S4) | no |
| 6 | `make_all_figures.py` — composite labels, 600 dpi, pixel-clamped | no |
| 7 | `build_tables.py` — Tables 1–7, Supplementary S1–S16 and Figures S1–S6, `RESULTS_NUMBERS.json` | no |
| 8 | `format_supp_tables.py` — publication formatting of S1–S16 | no |
| 9 | `presubmission_checks.py` | no |
| 10 | `verify_outputs.py` — the gate; the job exits non-zero if it fails | no |

No manuscript prose is generated.

## 5. Check before you trust it

```bash
L=$(ls -1dt outputs/pan_proteome_build_* | head -1)

python verify_outputs.py --inputs "$L"                 # expect 0 failures
cat "$L/tree_rooting.csv"                              # both families outgroup-rooted
ls "$L/final/supplementary/tables_formatted" | wc -l   # 16
head -3 "$L/final/supplementary/tables_formatted/TableS14_domain_classification.csv"
python -c "import json;N=json.load(open('$L/final/manuscript/RESULTS_NUMBERS.json'));print(len(N),'keys |',N['domains']['confidence'])"
```

Expect 24 JSON keys (including `rca_isoforms`) and `{'high confidence': 187}`. If `tree_rooting.csv` reports
midpoint for either family, the outgroup cache was cold — re-run `setup_env.sh` on a
login node and resubmit; the manuscript claims outgroup rooting in three places.

## 6. Send the results

```bash
bash run_all.sh collect
```

Bundles `final/`, `audit/`, `domains/`, `transcripts/`, `checks/`, `tree_rooting.csv`,
the core CSVs and `VERIFICATION.txt` with a manifest of sizes, timestamps and MD5s, and
prints the `scp` line.

---

## Table and figure numbering (manuscript v2.1)

Fixed by `numbering.py`, in order of first citation in the text:

| Main | Content | | Supplementary | File |
|---|---|---|---|---|
| Table 1 | summary, recruitment specificity | | S1 | `diamond_hits_all.csv` |
| Table 2 | protein models and gene loci per cultivar | | S2 | `audit/cpn60_classes.csv` |
| Table 3 | Chinese Spring gene loci (`cs_locus_table.csv`) | | S3 | `proteotypes.csv` |
| Table 4 | prevalence partition at three granularities | | S4 | `pav_matrix.csv` |
| Table 5 | **RCA proteotypes by established isoform (new)** | | S5 | `audit/novel_divergence.csv` |
| Table 6 | physicochemical descriptors | | S6 | `members.csv` |
| Table 7 | prioritised candidates (+ established isoform column) | | **S7** | **`audit/rca_isoform_assignment.csv` (new)** |
| | | | S8 | `proteotype_annotation.csv` |
| Figure S1 | RCA presence-absence matrix | | S9 | `polymorphic_sites.csv` |
| Figure S2 | Cpn60 presence-absence matrix | | S10 | `domains/domain_classification.csv` |
| **Figure S3** | **RCA residue-level alignment (new, PDF)** | | S11 | `transcripts/proteotype_recovery.csv` |
| **Figure S4** | **Cpn60 residue-level alignment (new, PDF)** | | S12 | `physicochemical_features.csv` |
| Figure S5 | genome-wide prevalence spectrum | | S13 | `proteotype_naming.csv` |
| Figure S6 | genome-wide per-cultivar composition | | S14 | `genome_wide_pav.csv` |
| | | | S15 | `pangenome_partition.csv` |
| | | | S16 | `group_comparison.csv` (empty: no phenotype labels) |

Each formatted file is:

```
"Bapat et al., 2026. <manuscript title>"
"Supplemental Table Sn. <table title>"
<column header row>
<data>
<blank row>
<column>,<description>          one row per column
```

`format_supp_tables.py` **fails** if a file contains a column with no description, or
carries a description for a column that is not in the file. Cross-references inside the
descriptions are resolved from `numbering.py`, so a renumbering cannot leave one stale.

## Fixes carried in this bundle

1. **Figure resolution** — matplotlib DPI set to 600 with a clamp keeping any figure
   under 60 megapixels and never below 300 dpi. A 141-row heatmap previously reached 189
   megapixels, which Pillow refuses to open.
2. **Panel letters** — the pipeline's pan-genome figure now draws (a)/(b), which its
   caption cites.
3. **Deterministic composite identifiers** — ties at the best e-value are broken on
   `(evalue, cultivar, source_id)`. Reordering the input previously changed 13 of 175
   names. Applied to both copies of the naming function.
4. **Silent midpoint rooting** — a failed outgroup lookup no longer falls back quietly;
   the run logs an error and writes `tree_rooting.csv` recording requested versus
   achieved rooting per family.
5. **Blank Table 1 cells** — `reference_query_sets.csv` lives outside the build
   directory; the lookup now falls back to a recursive search and warns if still empty.
6. **Re-runnable table build** — sources inside the output bundle are de-prioritised but
   not hidden, paths are compared by `realpath` so the `outputs/latest` symlink cannot
   cause a file to be copied onto itself, and destinations are made writable.
7. **Complete `RESULTS_NUMBERS.json`** — `domains` and `transcripts` sections are
   collected; they were absent before.
8. **Supplementary S14 and S15** — added to `SUPP_TABLES`; the manuscript cites them
   seven times and they were not being generated.
9. **Offline transcript step** — skips the Zenodo API entirely when the cache is complete.
10. **Relative project paths** — a run launched from this directory writes here, not into
    the parent.

## Still yours

Phenotype labels are unset under `cultivars:` in `config.yaml`, so Supplementary Table S9
is generated empty and the association framework is reported as a future application.
Setting two or more cultivars per group and re-running populates it with no other change.
