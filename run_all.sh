#!/usr/bin/env bash
# Orchestrator. Every step reports what it is doing and why it stopped.
#
#   bash run_all.sh setup     # LOGIN NODE, once: miniforge + env + data + caches
#   bash run_all.sh submit    # sbatch the job (does everything else)
#   bash run_all.sh status    # queue + tail of the log
#   bash run_all.sh finish [dir]  # run post-processing on an existing build
#   bash run_all.sh verify    # re-run the gate on the newest build
#   bash run_all.sh collect   # bundle the outputs for review
set -uo pipefail            # NOT -e: a failing step must report, not vanish

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
CONDA="$ROOT/miniforge3/etc/profile.d/conda.sh"
ENVD="$ROOT/envs/wheat_panproteome"

die () { echo; echo "ERROR: $*" >&2; exit 1; }

act () {
  if [[ ! -f "$CONDA" ]]; then
    if command -v conda >/dev/null 2>&1; then
      # shellcheck disable=SC1090
      source "$(conda info --base)/etc/profile.d/conda.sh"
    else
      die "conda not found at $CONDA and not on PATH. Run:  bash run_all.sh setup"
    fi
  else
    # shellcheck disable=SC1090
    source "$CONDA"
  fi
  [[ -x "$ENVD/bin/python" ]] || die "environment missing at $ENVD. Run:  bash run_all.sh setup"
  conda activate "$ENVD" || die "could not activate $ENVD"
  echo "[env] $(python -V 2>&1) from $ENVD"
}

latest () {
  local od l
  od="$(python - <<'PY' 2>/dev/null || echo outputs
import os, yaml
c = yaml.safe_load(open("config.yaml")) or {}
print(os.path.abspath((c.get("project") or {}).get("outputs_dir", "outputs")))
PY
)"
  l="$(ls -1dt "$od"/pan_proteome_build_* 2>/dev/null | head -1)"
  [[ -n "$l" ]] || die "no build directory under $od. Has the job run?"
  echo "$l"
}

case "${1:-help}" in

setup)
  echo ">>> LOGIN-NODE SETUP: miniforge, environment, proteomes, caches"
  bash setup_env.sh || die "setup_env.sh failed - see the output above"
  echo; echo ">>> setup complete. Next:  bash run_all.sh submit"
  ;;

submit)
  [[ -x "$ENVD/bin/python" ]] || die "no environment yet. Run:  bash run_all.sh setup"
  for d in inputs/hmm_cache inputs/zenodo_cache; do
    [[ -d "$d" ]] || echo "WARNING: $d missing; run 'act && python prefetch_caches.py' first"
  done
  command -v sbatch >/dev/null 2>&1 || die "sbatch not found - are you on a login node?"
  JID=$(sbatch --parsable submit_panproteome.sh) || die "sbatch rejected the job"
  echo ">>> submitted job $JID"
  echo "    watch:  tail -f logs/slurm_panproteome_${JID}.out"
  echo "    queue:  squeue -u \$USER"
  ;;

status)
  squeue -u "$USER" 2>/dev/null || echo "(squeue unavailable)"
  L=$(ls -1t logs/slurm_panproteome_*.out 2>/dev/null | head -1)
  [[ -n "${L:-}" ]] && { echo; echo "--- tail $L ---"; tail -30 "$L"; }
  ;;

finish)  shift || true; bash finish_build.sh "${1:-}" ;;

verify)  act; python verify_outputs.py --inputs "$(latest)" ;;

collect)
  act; L=$(latest); S="$ROOT/review_$(date +%Y%m%d_%H%M)"
  mkdir -p "$S"
  for d in final audit domains transcripts checks; do
    [[ -d "$L/$d" ]] && cp -r "$L/$d" "$S/"
  done
  cp "$L/tree_rooting.csv" "$S/" 2>/dev/null || true
  for f in proteotypes.csv members.csv pav_matrix.csv proteotype_naming.csv; do
    [[ -f "$L/$f" ]] && cp "$L/$f" "$S/"
  done
  python verify_outputs.py --inputs "$L" > "$S/VERIFICATION.txt" 2>&1
  tar czf "$S.tar.gz" -C "$ROOT" "$(basename "$S")" && rm -rf "$S"
  echo ">>> $S.tar.gz  ($(du -h "$S.tar.gz" | cut -f1))"
  echo "    scp kbapat@setonix.pawsey.org.au:$S.tar.gz ~/Downloads/"
  ;;

help|*)
  sed -n '2,9p' "$0"
  ;;
esac
