#!/bin/bash
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

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
cd "$PROJECT_DIR"
mkdir -p logs

echo "=================================================================="
echo " Wheat pan-proteome pipeline"
echo " Date        : $(date)"
echo " Host        : $(hostname)"
echo " Job ID      : ${SLURM_JOB_ID:-interactive}"
echo " Project dir : $PROJECT_DIR"
echo " CPUs        : ${SLURM_CPUS_PER_TASK:-?}"
echo "=================================================================="

for f in hpc_pipeline.py config.yaml environment.yml tree_tools.py external_data.py \
         panproteome_audit.py domain_profiling.py transcript_crossref.py rca_isoforms.py \
         alignment_views.py make_nature_figures.py build_tables.py; do
  [[ -f "$f" ]] || { echo "CRITICAL: $f not found in $PROJECT_DIR"; exit 1; }
done

for d in inputs/hmm_cache inputs/zenodo_cache; do
  [[ -d "$PROJECT_DIR/$d" ]] || echo "WARNING: $d absent - run 'python prefetch_caches.py' \
with internet access or the domain/transcript steps will fail."
done

ENV_PREFIX="${CONDA_ENV_PREFIX:-$PROJECT_DIR/envs/wheat_panproteome}"
CONDA_ROOT="${CONDA_ROOT:-$PROJECT_DIR/miniforge3}"

if [[ ! -x "$ENV_PREFIX/bin/python" ]]; then
  echo "CRITICAL: conda environment not found at $ENV_PREFIX"
  echo "  Run first:   bash setup_env.sh"
  exit 1
fi

CONDA_CACHE_ROOT="${CONDA_CACHE_ROOT:-$PROJECT_DIR/.conda_cache}"
export CONDA_PKGS_DIRS="$CONDA_CACHE_ROOT/pkgs"
export PIP_CACHE_DIR="$CONDA_CACHE_ROOT/pip"
export XDG_CACHE_HOME="$CONDA_CACHE_ROOT/xdg"
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$XDG_CACHE_HOME"

if [[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  source "$CONDA_ROOT/etc/profile.d/conda.sh"
  conda activate "$ENV_PREFIX"
else
  export PATH="$ENV_PREFIX/bin:$PATH"
fi
echo "[env] python : $(python -V 2>&1)  ($(command -v python))"
echo "[env] diamond: $(command -v diamond || echo MISSING)"
echo "[env] mafft  : $(command -v mafft   || echo MISSING)"
echo "[env] fasttree: $(command -v FastTree || command -v fasttree || echo MISSING)"

export TMPDIR="${TMPDIR:-/tmp}/wheat_pp_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
export MPLBACKEND=Agg
export MPLCONFIGDIR="$TMPDIR/mpl"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS="$THREADS"

python - <<'PY'
import os, yaml
c = yaml.safe_load(open("config.yaml")) or {}
p = c.get("project") or {}
i = os.path.abspath(p.get("inputs_dir", "inputs"))
o = os.path.abspath(p.get("outputs_dir", "outputs"))
print(f"[paths] launch dir : {os.path.abspath('.')}")
print(f"[paths] inputs_dir : {i}")
print(f"[paths] outputs_dir: {o}")
PY
python hpc_pipeline.py --check-env --config config.yaml || true

echo "[run] starting pipeline with $THREADS threads ..."
python hpc_pipeline.py --config config.yaml --threads "$THREADS"

OUTDIR="$(python - <<'PY'
import os, yaml
c = yaml.safe_load(open("config.yaml")) or {}
print(os.path.abspath((c.get("project") or {}).get("outputs_dir", "outputs")))
PY
)"
LATEST="$(ls -1dt "$OUTDIR"/pan_proteome_build_* 2>/dev/null | head -1 || true)"
if [[ -z "${LATEST:-}" ]]; then
  echo "CRITICAL: no result directory under $OUTDIR"
  exit 1
fi
echo "[run] result: $LATEST"
ln -sfn "$LATEST" "$OUTDIR/latest"

echo "[post] 1/8 locus resolution, specificity, classes, divergence, granularity ..."
python panproteome_audit.py --inputs "$LATEST" --outdir "$LATEST/audit" || true

echo "[post] 2/8 domain profiling (HMM + reciprocal best hits) ..."
python domain_profiling.py --inputs "$LATEST" --outdir "$LATEST/domains" \
    --cache "$PROJECT_DIR/inputs/hmm_cache" --threads "$THREADS" \
  || echo "[post] WARNING domain profiling failed - was prefetch_caches.py run?"

echo "[post] 3/8 cross-reference to the independent annotation ..."
python transcript_crossref.py --inputs "$LATEST" --outdir "$LATEST/transcripts" \
    --cache "$PROJECT_DIR/inputs/zenodo_cache" \
  || echo "[post] WARNING transcript cross-reference failed - was prefetch_caches.py run?"

echo "[post] 4/8 locus resolution again, after domain profiling ..."
python panproteome_audit.py --inputs "$LATEST" --outdir "$LATEST/audit" || true

echo "[post] 5/8 RCA isoform assignment ..."
python rca_isoforms.py --inputs "$LATEST" --outdir "$LATEST/audit" --config config.yaml \
  || echo "[post] WARNING rca_isoforms.py failed"

echo "[post] 6/8 residue-level alignment views ..."
python alignment_views.py --inputs "$LATEST" --outdir "$LATEST/figures" \
  || echo "[post] WARNING alignment_views.py failed"

echo "[post] 7/8 figures and source data ..."
python make_nature_figures.py --inputs "$LATEST" --outdir "$LATEST/figures/nature" \
    --source-data "$LATEST/source_data" \
  || echo "[post] WARNING make_nature_figures.py failed"

echo "[post] 8/8 tables, supplementary data and output bundle ..."
python build_tables.py --inputs "$LATEST" --outdir "$LATEST/final" --config config.yaml \
  || echo "[post] WARNING build_tables.py failed"

echo "=================================================================="
echo " Completed : $(date)"
echo " Build     : $LATEST"
echo " Outputs   : $LATEST/final"
echo "=================================================================="
