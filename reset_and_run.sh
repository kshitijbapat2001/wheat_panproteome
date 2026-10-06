#!/usr/bin/env bash
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPECTED="/scratch/pawsey1276/kbapat/Pan_prot_final/FINAL"
STAMP="$(date +%Y%m%d_%H%M)"
PLAN=0; PURGE_INPUTS=0; ASSUME_YES=0
for a in "$@"; do
  case "$a" in
    --plan) PLAN=1 ;;
    --purge-inputs) PURGE_INPUTS=1 ;;
    --yes|-y) ASSUME_YES=1 ;;
    *) echo "unknown option: $a"; exit 2 ;;
  esac
done
cd "$ROOT" || exit 1

die () { echo; echo "ERROR: $*" >&2; exit 1; }
hr  () { printf '%s\n' "----------------------------------------------------------------------"; }

if [[ "$ROOT" != "$EXPECTED" ]]; then
  echo "WARNING: running in $ROOT, not the expected $EXPECTED"
  if [[ $PLAN -eq 0 && $ASSUME_YES -eq 0 ]]; then
    if [[ -t 0 ]]; then
      read -rp "continue anyway? [y/N] " a; [[ "$a" == [yY] ]] || exit 1
    else
      die "not an interactive shell; re-run with --yes if this directory is correct"
    fi
  fi
fi
[[ -f "$ROOT/hpc_pipeline.py" && -f "$ROOT/config.yaml" ]] \
  || die "this does not look like the pipeline directory (hpc_pipeline.py / config.yaml missing)"

REMOVE=(envs miniforge3 outputs logs __pycache__ github_stage zenodo_upload
        supplementary_tables_formatted .conda_cache)
[[ $PURGE_INPUTS -eq 1 ]] && REMOVE+=(inputs)

hr; echo " RESET PLAN"; echo "   directory : $ROOT"; hr
echo " will REMOVE:"
for d in "${REMOVE[@]}"; do
  [[ -e "$ROOT/$d" ]] && printf "   %-32s %s\n" "$d/" "$(du -sh "$ROOT/$d" 2>/dev/null | cut -f1)"
done
for g in "review_*" "*.tar.gz" "*.zip" "*.bak" "*.bak_ids" "*.pre_v3" "*.orig.png"; do
  ls -d $ROOT/$g >/dev/null 2>&1 && printf "   %-32s\n" "$g"
done
echo
echo " will KEEP:"
printf "   %-32s %s\n" "code (*.py *.sh *.yaml *.yml *.md)" "$(ls "$ROOT"/*.py "$ROOT"/*.sh "$ROOT"/*.yaml "$ROOT"/*.yml "$ROOT"/*.md 2>/dev/null | wc -l) files"
if [[ $PURGE_INPUTS -eq 0 && -d "$ROOT/inputs" ]]; then
  printf "   %-32s %s  (proteomes and warmed caches; --purge-inputs to delete)\n" \
         "inputs/" "$(du -sh "$ROOT/inputs" 2>/dev/null | cut -f1)"
fi
[[ -d "$ROOT/outputs" ]] && echo && echo " existing outputs/ will be ARCHIVED to archive/outputs_$STAMP.tar.gz first"
hr

[[ $PLAN -eq 1 ]] && { echo "--plan: nothing changed."; exit 0; }
if [[ $ASSUME_YES -eq 0 ]]; then
  [[ -t 0 ]] || die "not an interactive shell; re-run with --yes to confirm the reset"
  read -rp "Type RESET to proceed: " ans
  [[ "$ans" == "RESET" ]] || { echo "aborted."; exit 1; }
fi

if [[ -d "$ROOT/outputs" ]]; then
  mkdir -p "$ROOT/archive"
  echo "[reset] archiving existing outputs ..."
  tar czf "$ROOT/archive/outputs_$STAMP.tar.gz" -C "$ROOT" outputs \
    && echo "[reset] archived -> archive/outputs_$STAMP.tar.gz ($(du -h "$ROOT/archive/outputs_$STAMP.tar.gz" | cut -f1))" \
    || die "archiving failed - refusing to delete outputs/"
fi
for d in "${REMOVE[@]}"; do
  [[ -n "$d" && -e "$ROOT/$d" ]] && { echo "[reset] removing $d/"; rm -rf "${ROOT:?}/$d"; }
done
for g in review_ ; do rm -rf "${ROOT:?}"/${g}* 2>/dev/null; done
rm -f "${ROOT:?}"/*.tar.gz "${ROOT:?}"/*.zip 2>/dev/null
rm -f "${ROOT:?}"/*.bak "${ROOT:?}"/*.bak_ids "${ROOT:?}"/*.pre_v3 2>/dev/null
find "${ROOT:?}" -name "*.orig.png" -delete 2>/dev/null
find "${ROOT:?}" -name "__pycache__" -type d -prune -exec rm -rf {} + 2>/dev/null
mkdir -p "$ROOT/logs"
echo "[reset] done."; hr

echo "[1/4] environment, proteomes, references and caches (login node, needs internet)"
bash setup_env.sh 2>&1 | tee "logs/setup_$STAMP.log" | tail -20
[[ -x "$ROOT/envs/wheat_panproteome/bin/python" ]] || die "environment was not created - see logs/setup_$STAMP.log"

source "$ROOT/miniforge3/etc/profile.d/conda.sh" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ROOT/envs/wheat_panproteome" || die "cannot activate the environment"

for d in inputs/hmm_cache inputs/zenodo_cache; do
  [[ -d "$ROOT/$d" ]] || die "$d missing after setup - the job cannot run the domain/transcript steps offline"
done
echo "[1/4] OK"; hr

echo "[2/4] submitting the pipeline job"
command -v sbatch >/dev/null 2>&1 || die "sbatch not found - are you on a login node?"
JID=$(sbatch --parsable submit_panproteome.sh) || die "sbatch rejected the job"
echo "      job $JID submitted; waiting (tail -f logs/slurm_panproteome_${JID}.out)"
while squeue -j "$JID" -h 2>/dev/null | grep -q .; do sleep 60; done
echo "[2/4] job finished"; hr

echo "[3/4] verification"
L="$(ls -1dt "$ROOT"/outputs/pan_proteome_build_* 2>/dev/null | head -1)"
[[ -n "$L" ]] || die "no build directory produced - see logs/slurm_panproteome_${JID}.err"
python verify_outputs.py --inputs "$L"; RC=$?
hr

echo "[4/4] publication-formatted supplementary tables"
python format_supp_tables.py "$L/final/supplementary/tables" "$L/final/supplementary/tables_formatted"
echo
echo "======================================================================"
echo " build      : $L"
echo " outputs    : $L/final"
echo " supp tables: $L/final/supplementary/tables_formatted   (S1-S16, publication format)"
echo " archive    : archive/outputs_$STAMP.tar.gz"
if [[ $RC -ne 0 ]]; then
  echo " *** VERIFICATION FAILED - do not use these outputs until resolved ***"
else
  echo " Verification passed. Bundle for upload with:  bash run_all.sh collect"
fi
echo "======================================================================"
exit $RC
