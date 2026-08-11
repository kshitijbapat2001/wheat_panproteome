# Wheat RCA / Cpn60 pan-proteome pipeline

Builds a gene-locus-anchored catalogue of Rubisco activase and chaperonin-60 protein
variation across ten sequenced bread wheat cultivars, from Ensembl Plants proteomes.

## Run

```bash
bash run_all.sh setup     # login node: miniforge, environment, proteomes, caches
bash run_all.sh submit    # SLURM job: everything else
bash run_all.sh verify    # consistency gate
```

`RUNBOOK.md` has the detail, including the login-node/compute-node split — compute nodes
have no internet, so all downloads are cached beforehand by `prefetch_caches.py`.

## What it produces

Proteotype catalogue with composite identifiers, gene-locus resolution, presence–absence
matrices, recruitment-specificity audit, Cpn60 subunit classes, HMM domain profiling with
an out-of-family negative control, cross-reference against an independent de novo
annotation, alignments and rooted phylogenies, genome-wide pan-proteome context, and the
main and supplementary tables and figures. No manuscript prose is generated.

## Verification

`verify_outputs.py` gates a build and fails on any of: mixed builds, disagreeing
proteotype counts, outgroup rooting requested but not achieved, a truncated candidate
table, blank cells in Table 1, figures below 300 dpi or above 60 megapixels, or malformed
composite identifiers. A non-zero exit means the outputs should not be quoted.

## Data

Derived datasets are deposited separately under restricted access; see the Data
availability statement in the associated manuscript.

## Licence

MIT (code). Derived data are released under CC BY 4.0.
