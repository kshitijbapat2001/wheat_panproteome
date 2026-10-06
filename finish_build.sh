#!/usr/bin/env bash
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
CONDA="$ROOT/miniforge3/etc/profile.d/conda.sh"
ENVD="$ROOT/envs/wheat_panproteome"
die () { echo; echo "ERROR: $*" >&2; exit 1; }

if [[ -f "$CONDA" ]]; then source "$CONDA"; elif command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"; else die "conda not found"; fi
[[ -x "$ENVD/bin/python" ]] || die "environment missing at $ENVD"
conda activate "$ENVD" || die "cannot activate $ENVD"

OUTDIR="$(python - <<'PY'
import os, yaml
c = yaml.safe_load(open("config.yaml")) or {}
print(os.path.abspath((c.get("project") or {}).get("outputs_dir", "outputs")))
PY
)"
L="${1:-}"
if [[ -z "$L" ]]; then
  L="$(ls -1dt "$OUTDIR"/pan_proteome_build_* 2>/dev/null | head -1)"
  [[ -n "$L" ]] || L="$(ls -1dt "$ROOT"/../*/outputs/pan_proteome_build_* 2>/dev/null | head -1)"
fi
[[ -n "$L" && -d "$L" ]] || die "no build directory found. Pass one explicitly:
    bash finish_build.sh /scratch/.../outputs/pan_proteome_build_YYYYMMDD_HHMMSS"
L="$(readlink -f "$L")"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "=================================================================="
echo " Finishing build"
echo "   build   : $L"
echo "   threads : $THREADS"
echo "=================================================================="
[[ -f "$L/proteotypes.csv" ]] || die "$L does not look like a completed build (no proteotypes.csv)"

echo "[post] 1/6 audit ..."
python panproteome_audit.py --inputs "$L" --outdir "$L/audit" || true
echo "[post] 2/6 domain profiling ..."
python domain_profiling.py --inputs "$L" --outdir "$L/domains" \
    --cache "$ROOT/inputs/hmm_cache" --threads "$THREADS" \
  || echo "[post] WARNING domain profiling failed (needs internet or a warmed hmm_cache)"
echo "[post] 3/6 transcript cross-reference ..."
python transcript_crossref.py --inputs "$L" --outdir "$L/transcripts" \
    --cache "$ROOT/inputs/zenodo_cache" \
  || echo "[post] WARNING transcript cross-reference failed (needs internet or a warmed cache)"
echo "[post] 4/6 audit again (domain classes now exist) ..."
python panproteome_audit.py --inputs "$L" --outdir "$L/audit" || true
echo "[post] 4b RCA isoform assignment (rca_isoforms.py) ..."
python rca_isoforms.py --inputs "$L" --outdir "$L/audit" --config config.yaml \
  || echo "[post] WARNING rca_isoforms.py failed - Table 5, Supplementary Table S7 and Figure S3 need it"
echo "[post] 4c residue-level alignment views (alignment_views.py) ..."
python alignment_views.py --inputs "$L" --outdir "$L/figures" \
  || echo "[post] WARNING alignment_views.py failed - Supplementary Figures S3/S4 need it"
echo "[post] 5/6 figures ..."
python make_all_figures.py --inputs "$L" --outdir "$L/final/figures" \
    --suppdir "$L/final/supplementary" --label composite --cpn60-rows 30 \
    --sequences "$L/PanProteome_catalog.faa" \
  || python make_all_figures.py --inputs "$L" --outdir "$L/final/figures" \
        --suppdir "$L/final/supplementary" --label composite --cpn60-rows 30 || true
echo "[post] 5b Nature-format figures and source data (make_nature_figures.py) ..."
python make_nature_figures.py --inputs "$L" --outdir "$L/final/nature_figures" \
    --source-data "$L/final/source_data" \
  || echo "[post] WARNING make_nature_figures.py failed"
echo "[post] 6/6 tables + RESULTS_NUMBERS.json ..."
python build_tables.py --inputs "$L" --outdir "$L/final" --config config.yaml || true
python presubmission_checks.py --inputs "$L" --outdir "$L/checks" || true
echo "[post] publication-formatted supplementary tables ..."
python format_supp_tables.py "$L/final/supplementary/tables" \
                             "$L/final/supplementary/tables_formatted" || true

mkdir -p "$OUTDIR"; ln -sfn "$L" "$OUTDIR/latest"
echo "=================================================================="
python verify_outputs.py --inputs "$L"; RC=$?
echo "=================================================================="
echo " build   : $L"
echo " outputs : $L/final"
[[ $RC -eq 0 ]] && echo " Verification passed." || echo " *** VERIFICATION FAILED - fix before use. ***"
exit $RC
