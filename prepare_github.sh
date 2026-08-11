#!/usr/bin/env bash
# Stage the pipeline code for the private GitHub repository, and build the code archive
# for Zenodo. Code only: no proteomes, no build outputs, no caches, no environments.
#
#   bash prepare_github.sh
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
STAGE="$ROOT/github_stage"; rm -rf "$STAGE"; mkdir -p "$STAGE"

CODE=(hpc_pipeline.py panproteome_audit.py tree_tools.py external_data.py
      domain_profiling.py transcript_crossref.py make_all_figures.py
      build_tables.py presubmission_checks.py verify_outputs.py
      prefetch_caches.py diagnose_inputs.py shrink_oversized_figures.py
      rebuild_naming.py regen_id_figures.py regen_lowres_figures.py
      check_figure_dpi.py apply_v3_fixes.py whats_fresh.sh
      run_all.sh setup_env.sh stage_data.sh submit_panproteome.sh finish_build.sh
      prepare_zenodo.sh prepare_github.sh
      config.yaml environment.yml RUNBOOK.md)

n=0
for f in "${CODE[@]}"; do
  if [[ -f "$f" ]]; then cp "$f" "$STAGE/"; n=$((n+1)); else echo "  - $f (absent)"; fi
done
echo "  staged $n file(s)"

cat > "$STAGE/.gitignore" <<'GI'
# data and results never belong in the repository
outputs/
inputs/
envs/
miniforge3/
logs/
zenodo_upload/
github_stage/
review_*/
*.tar.gz
*.zip
# large or generated artefacts
*.faa
*.fasta
*.nwk
*.aln
*.dmnd
*.hmm
*.db
*.sqlite
*.png
*.pdf
*.docx
# backups the patch scripts leave behind
*.bak
*.bak_ids
*.pre_v3
*.orig.png
__pycache__/
*.pyc
.conda_cache/
GI

cat > "$STAGE/README.md" <<'MD'
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
MD

cat > "$STAGE/LICENSE" <<'LC'
MIT License

Copyright (c) 2026 Kshitij Bapat, Nicolas Taylor, Agyeya Pratap

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
LC

TAR="$ROOT/wheat-panproteome-code-$(date +%Y%m%d).tar.gz"
tar czf "$TAR" -C "$STAGE" .
echo
echo "  staged repo : $STAGE"
echo "  code archive: $TAR  ($(du -h "$TAR" | cut -f1))   <- upload this to Zenodo"
echo
echo "Sanity check - nothing large or private should appear below:"
find "$STAGE" -type f -size +1M -printf "  LARGE %s %p\n" 2>/dev/null || true
echo "  total $(find "$STAGE" -type f | wc -l) files, $(du -sh "$STAGE" | cut -f1)"
