#!/bin/bash
# =====================================================================================
# submit_panproteome.sh  --  SLURM launcher, Pawsey Setonix
#
#   cd /scratch/pawsey1276/kbapat/Pan_prot_final
#   bash setup_env.sh          # ONCE, on a login node (internet + conda env + data)
#   sbatch submit_panproteome.sh
#
# Setonix compute nodes have no outbound internet, so setup_env.sh must run first on a
# login node: it builds the conda environment and pre-downloads every proteome, every
# UniProt reference and every external annotation into inputs/. This job then runs
# entirely offline.
# =====================================================================================
#SBATCH --job-name=wheat_panproteome
#SBATCH --account=pawsey1276
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=115G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm_panproteome_%j.out
#SBATCH --error=logs/slurm_panproteome_%j.err
#SBATCH --export=NONE
# Resource notes for Setonix `work`: a node is 128 cores / ~230 GB, allocated
# proportionally. 64 cores -> ~115 GB, which is the half-node request above. The
# whole-proteome stage is the expensive one; DIAMOND memory is roughly
# 6 x whole_proteome.block_size_gb, and whole_proteome.memory_limit_gb must stay well
# below --mem. A DIAMOND step dying with rc=-9 was killed by the OOM killer: lower
# block_size_gb in config.yaml, or request a full node (--exclusive --cpus-per-task=128).
# The stage caches its per-cultivar searches (whole_proteome.resume: true), so a job that
# hits the wall clock can simply be re-submitted and will resume.

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-/scratch/pawsey1276/kbapat/Pan_prot_final}}"
cd "$PROJECT_DIR"
mkdir -p logs

echo "=================================================================="
echo " Wheat pan-proteome pipeline (Setonix)"
echo " Date        : $(date)"
echo " Host        : $(hostname)"
echo " Job ID      : ${SLURM_JOB_ID:-interactive}"
echo " Project dir : $PROJECT_DIR"
echo " CPUs        : ${SLURM_CPUS_PER_TASK:-?}"
echo "=================================================================="

for f in hpc_pipeline.py config.yaml environment.yml tree_tools.py \
         panproteome_audit.py external_data.py build_tables.py \
         make_all_figures.py domain_profiling.py transcript_crossref.py \
         presubmission_checks.py verify_outputs.py; do
  [[ -f "$f" ]] || { echo "CRITICAL: $f not found in $PROJECT_DIR"; exit 1; }
done

# Caches must have been warmed on a login node: compute nodes have no internet.
for d in inputs/hmm_cache inputs/zenodo_cache; do
  [[ -d "$PROJECT_DIR/$d" ]] || echo "WARNING: $d absent - run 'python prefetch_caches.py' \
on a login node or the domain/transcript steps will fail."
done

# --- environment ------------------------------------------------------------------
ENV_PREFIX="${CONDA_ENV_PREFIX:-$PROJECT_DIR/envs/wheat_panproteome}"
CONDA_ROOT="${CONDA_ROOT:-$PROJECT_DIR/miniforge3}"

if [[ ! -x "$ENV_PREFIX/bin/python" ]]; then
  echo "CRITICAL: conda environment not found at $ENV_PREFIX"
  echo "  Run this on a LOGIN node first:   bash setup_env.sh"
  exit 1
fi

# Caches off HOME (small quota on Setonix); everything under the project dir.
CONDA_CACHE_ROOT="${CONDA_CACHE_ROOT:-$PROJECT_DIR/.conda_cache}"
export CONDA_PKGS_DIRS="$CONDA_CACHE_ROOT/pkgs"
export PIP_CACHE_DIR="$CONDA_CACHE_ROOT/pip"
export XDG_CACHE_HOME="$CONDA_CACHE_ROOT/xdg"
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$XDG_CACHE_HOME"

if [[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  # shellcheck disable=SC1091
  source "$CONDA_ROOT/etc/profile.d/conda.sh"
  conda activate "$ENV_PREFIX"
else
  # Fall back to using the env's bin directly (works without `conda activate`).
  export PATH="$ENV_PREFIX/bin:$PATH"
fi
echo "[env] python : $(python -V 2>&1)  ($(command -v python))"
echo "[env] diamond: $(command -v diamond || echo MISSING)"
echo "[env] mafft  : $(command -v mafft   || echo MISSING)"
echo "[env] fasttree: $(command -v FastTree || command -v fasttree || echo MISSING)"

# Node-local scratch for tool temp files; all outputs stay in the project dir.
export TMPDIR="${TMPDIR:-/tmp}/wheat_pp_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
export MPLBACKEND=Agg
export MPLCONFIGDIR="$TMPDIR/mpl"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS="$THREADS"

# --- pre-flight --------------------------------------------------------------------
python - <<'PY'
import os, yaml
c = yaml.safe_load(open("config.yaml")) or {}
p = c.get("project") or {}
i = os.path.abspath(p.get("inputs_dir", "inputs"))
o = os.path.abspath(p.get("outputs_dir", "outputs"))
here = os.path.abspath(".")
print(f"[paths] launch dir : {here}")
print(f"[paths] inputs_dir : {i}")
print(f"[paths] outputs_dir: {o}")
for lbl, d in (("inputs", i), ("outputs", o)):
    if not d.startswith(here):
        print(f"[paths] NOTE {lbl}_dir is OUTSIDE the launch directory. That is allowed "
              f"(useful for sharing downloads) but results will not appear under {here}.")
PY
python hpc_pipeline.py --check-env --config config.yaml || true

# --- run ---------------------------------------------------------------------------
echo "[run] starting pipeline with $THREADS threads ..."
python hpc_pipeline.py --config config.yaml --threads "$THREADS"

# Resolve the results directory from config.yaml rather than assuming ./outputs. A
# hard-coded absolute outputs_dir once sent a run launched from FINAL/ into the parent
# project's outputs/, and the job then declared "no result directory produced" on a run
# that had in fact succeeded.
OUTDIR="$(python - <<'PY'
import os, yaml
c = yaml.safe_load(open("config.yaml")) or {}
print(os.path.abspath((c.get("project") or {}).get("outputs_dir", "outputs")))
PY
)"
LATEST="$(ls -1dt "$OUTDIR"/pan_proteome_build_* 2>/dev/null | head -1 || true)"
if [[ -z "${LATEST:-}" ]]; then
  echo "CRITICAL: no result directory under $OUTDIR"
  echo "          config.yaml project.outputs_dir resolves to: $OUTDIR"
  echo "          if the build went elsewhere, finish it with:"
  echo "              bash finish_build.sh /path/to/pan_proteome_build_YYYYMMDD_HHMMSS"
  exit 1
fi
echo "[run] result: $LATEST"

# --- post-processing, all offline from the caches warmed by setup_env.sh ------------
ln -sfn "$LATEST" "$PROJECT_DIR/outputs/latest"

echo "[post] 1/6 audit ..."
python panproteome_audit.py --inputs "$LATEST" --outdir "$LATEST/audit" || true

echo "[post] 2/6 domain profiling (HMM + reciprocal best hits) ..."
python domain_profiling.py --inputs "$LATEST" --outdir "$LATEST/domains" \
    --cache "$PROJECT_DIR/inputs/hmm_cache" --threads "$THREADS" \
  || echo "[post] WARNING domain profiling failed - was prefetch_caches.py run?"

echo "[post] 3/6 transcript cross-reference ..."
python transcript_crossref.py --inputs "$LATEST" --outdir "$LATEST/transcripts" \
    --cache "$PROJECT_DIR/inputs/zenodo_cache" \
  || echo "[post] WARNING transcript cross-reference failed - was prefetch_caches.py run?"

echo "[post] 4/6 audit again (now that domain classes exist) ..."
python panproteome_audit.py --inputs "$LATEST" --outdir "$LATEST/audit" || true

echo "[post] 5/6 figures (composite labels, 600 dpi) ..."
python make_all_figures.py --inputs "$LATEST" --outdir "$LATEST/final/figures" \
    --suppdir "$LATEST/final/supplementary" --label composite --cpn60-rows 30 \
    --sequences "$LATEST/PanProteome_catalog.faa" \
  || python make_all_figures.py --inputs "$LATEST" --outdir "$LATEST/final/figures" \
        --suppdir "$LATEST/final/supplementary" --label composite --cpn60-rows 30 || true

echo "[post] 6/6 tables + RESULTS_NUMBERS.json ..."
python build_tables.py --inputs "$LATEST" --outdir "$LATEST/final" \
                       --config config.yaml || true

echo "[post] presubmission checks ..."
python presubmission_checks.py --inputs "$LATEST" --outdir "$LATEST/checks" || true

# --- gate ---------------------------------------------------------------------------
echo "=================================================================="
echo " VERIFICATION"
echo "=================================================================="
python verify_outputs.py --inputs "$LATEST"
VERIFY_RC=$?

echo "=================================================================="
echo " Completed: $(date)"
echo " Result dir         : $LATEST"
echo " Stable symlink     : $PROJECT_DIR/outputs/latest"
echo " FINAL BUNDLE       : $LATEST/final/"
echo "   manuscript draft : $LATEST/final/manuscript/Manuscript_draft.docx"
echo "   main figures     : $LATEST/final/figures/"
echo "   main tables      : $LATEST/final/tables/"
echo "   supplementary    : $LATEST/final/supplementary/"
echo "   manifest         : $LATEST/final/MANIFEST.md"
echo " Audit tables       : $LATEST/audit/"
echo "=================================================================="
echo " Outputs  : $LATEST/final"
if [[ "${VERIFY_RC:-1}" -ne 0 ]]; then
  echo " *** VERIFICATION FAILED - see the block above before using these outputs. ***"
  exit "$VERIFY_RC"
fi
echo " Verification passed; outputs are safe to use."
