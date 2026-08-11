#!/bin/bash
# =====================================================================================
# setup_env.sh  --  ONE-TIME setup on a Setonix LOGIN node (needs internet)
#
#   cd /scratch/pawsey1276/kbapat/Pan_prot_final
#   bash setup_env.sh
#   sbatch submit_panproteome.sh
#
# Does four things, each idempotent, so it is safe to re-run:
#   1. installs Miniforge into the project directory if conda is not already available
#      (Setonix has no conda module, and HOME is too small for conda packages);
#   2. creates the pipeline conda environment from environment.yml;
#   3. pre-downloads every cultivar proteome and every UniProt reference into inputs/;
#   4. warms the external-database cache (UniProt/InterPro outgroups and annotation).
#
# After this, compute nodes never need the internet.
# =====================================================================================
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"
mkdir -p logs inputs outputs

echo "=================================================================="
echo " Setup: wheat pan-proteome pipeline"
echo " Project dir : $PROJECT_DIR"
echo " Host        : $(hostname)"
echo "=================================================================="

for f in hpc_pipeline.py config.yaml environment.yml tree_tools.py \
         panproteome_audit.py external_data.py build_tables.py \
         make_all_figures.py domain_profiling.py transcript_crossref.py \
         prefetch_caches.py verify_outputs.py; do
  [[ -f "$f" ]] || { echo "CRITICAL: $f missing from $PROJECT_DIR"; exit 1; }
done
# This build deliberately ships no manuscript generator; tables and numbers only.

# --- caches off HOME (Setonix HOME quota is small) ---------------------------------
CONDA_CACHE_ROOT="${CONDA_CACHE_ROOT:-$PROJECT_DIR/.conda_cache}"
export CONDA_PKGS_DIRS="$CONDA_CACHE_ROOT/pkgs"
export PIP_CACHE_DIR="$CONDA_CACHE_ROOT/pip"
export XDG_CACHE_HOME="$CONDA_CACHE_ROOT/xdg"
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$XDG_CACHE_HOME"

# --- 1. conda ----------------------------------------------------------------------
CONDA_ROOT="${CONDA_ROOT:-$PROJECT_DIR/miniforge3}"
if [[ ! -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  if command -v conda >/dev/null 2>&1; then
    CONDA_ROOT="$(conda info --base)"
    echo "[env] using conda already on PATH: $CONDA_ROOT"
  else
    echo "[env] installing Miniforge into $CONDA_ROOT (no conda module on Setonix) ..."
    INSTALLER="$PROJECT_DIR/.miniforge_installer.sh"
    URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    if command -v curl >/dev/null 2>&1; then
      curl -fsSL "$URL" -o "$INSTALLER"
    else
      wget -q "$URL" -O "$INSTALLER"
    fi
    bash "$INSTALLER" -b -p "$CONDA_ROOT"
    rm -f "$INSTALLER"
    echo "[env] Miniforge installed."
  fi
fi
# shellcheck disable=SC1091
source "$CONDA_ROOT/etc/profile.d/conda.sh"
SOLVER="conda"
command -v mamba >/dev/null 2>&1 && SOLVER="mamba"
echo "[env] solver: $SOLVER   base: $CONDA_ROOT"

# --- 2. environment ----------------------------------------------------------------
ENV_PREFIX="${CONDA_ENV_PREFIX:-$PROJECT_DIR/envs/wheat_panproteome}"
if [[ -x "$ENV_PREFIX/bin/python" ]]; then
  echo "[env] reusing existing environment: $ENV_PREFIX"
else
  echo "[env] creating environment: $ENV_PREFIX"
  "$SOLVER" env create --prefix "$ENV_PREFIX" --file environment.yml
fi
conda activate "$ENV_PREFIX"
echo "[env] python  : $(python -V 2>&1)"
for t in diamond mafft FastTree; do
  printf '[env] %-8s: %s\n' "$t" "$(command -v "$t" || echo MISSING)"
done

python hpc_pipeline.py --check-env --config config.yaml || true

# --- 3. proteomes + UniProt references ---------------------------------------------
echo "[stage] downloading proteomes and reference sequences (this can take a while) ..."
python - <<'PY'
import os, sys, yaml
sys.path.insert(0, os.getcwd())
import hpc_pipeline as H
cfg = H.load_config("config.yaml")
log = H.setup_logging(os.path.abspath(cfg["project"]["outputs_dir"]))
H._install_proxy(cfg, log)
inputs_dir = os.path.abspath(cfg["project"]["inputs_dir"])
os.makedirs(inputs_dir, exist_ok=True)
refs = H.fetch_reference_sequences(cfg, os.path.join(inputs_dir, "references"), log)
for fam in H.FAMILIES:
    n = sum(1 for f in refs["qid2family"].values() if f == fam)
    print(f"[stage] {fam}: {n} reference sequences")
base, rel = H.resolve_release(cfg, log)
prot = H.fetch_cultivar_proteomes(cfg, inputs_dir, base, rel, log)
print(f"[stage] proteomes staged: {len(prot)}/{len(cfg['cultivars'])} -> {inputs_dir}/proteomes")
PY

# --- 4. external database cache ----------------------------------------------------
echo "[stage] warming the UniProt/InterPro cache ..."
python external_data.py --warm-cache --config config.yaml || \
  echo "[stage] WARNING: external cache warm-up failed; trees will be midpoint-rooted."

# --- 5. HMM + Zenodo caches (compute nodes have no internet) ------------------------
echo "[stage] pre-fetching Pfam HMMs, Swiss-Prot seeds and the Zenodo tables ..."
python prefetch_caches.py --inputs "$PROJECT_DIR/inputs" || \
  echo "[stage] WARNING: prefetch incomplete; domain/transcript steps will fail in the job."

echo "=================================================================="
echo " Setup complete."
echo "   environment : $ENV_PREFIX"
echo "   inputs      : $PROJECT_DIR/inputs"
echo " Next:  sbatch submit_panproteome.sh"
echo "=================================================================="
