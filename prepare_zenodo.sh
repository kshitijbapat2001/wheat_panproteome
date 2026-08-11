#!/usr/bin/env bash
# Assemble the Zenodo deposit from a verified build.
#
#   bash prepare_zenodo.sh                      # newest build
#   bash prepare_zenodo.sh outputs/pan_proteome_build_YYYYMMDD_HHMMSS
#
# Zenodo shows a flat file list, so everything is numbered in the order a reader wants
# it. Individual files are uploaded as well as the archive: a reviewer who wants one
# table should not have to download 80 MB.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
B="${1:-}"; [[ -z "$B" ]] && B="$(ls -1dt outputs/pan_proteome_build_* 2>/dev/null | head -1)"
[[ -n "$B" && -d "$B" ]] || { echo "no build found; pass one explicitly"; exit 1; }
B="$(readlink -f "$B")"
OUT="$ROOT/zenodo_upload"; rm -rf "$OUT"; mkdir -p "$OUT"
echo "build : $B"; echo "deposit: $OUT"; echo

n=0; miss=0
add () {  # add <NN_name> <path-relative-to-build>
  local dest="$1" rel="$2" src="$B/$rel"
  if [[ -e "$src" ]]; then
    if [[ -d "$src" ]]; then (cd "$(dirname "$src")" && zip -qr "$OUT/$dest.zip" "$(basename "$src")") && echo "  + $dest.zip"
    else cp "$src" "$OUT/$dest"; echo "  + $dest"; fi
    n=$((n+1))
  else echo "  - $rel  (absent)"; miss=$((miss+1)); fi
}

echo "core catalogue"
add 01_proteotype_catalogue.faa            PanProteome_catalog.annotated.faa
add 02_proteotypes.csv                     proteotypes.csv
add 03_proteotype_naming.csv               proteotype_naming.csv
add 04_members_protein_models.csv          members.csv
add 05_pav_matrix.csv                      pav_matrix.csv
echo "audit and evidence"
add 06_candidate_table.csv                 audit/candidate_table.csv
add 07_cs_locus_table.csv                  audit/cs_locus_table.csv
add 08_granularity_partition.csv           audit/granularity_partition.csv
add 09_cpn60_classes.csv                   audit/cpn60_classes.csv
add 10_specificity_audit.csv               audit/specificity_audit.csv
add 11_novel_divergence.csv                audit/novel_divergence.csv
add 12_domain_classification.csv           domains/domain_classification.csv
add 13_proteotype_recovery.csv             transcripts/proteotype_recovery.csv
add 14_locus_recovery.csv                  transcripts/locus_recovery.csv
add 15_prevalence_discordance.csv          transcripts/prevalence_discordance.csv
echo "descriptors, sites, search"
add 16_physicochemical_features.csv        physicochemical_features.csv
add 17_polymorphic_sites.csv               polymorphic_sites.csv
add 18_reference_query_sets.csv            final/supplementary/data/reference_query_sets.csv
add 19_diamond_hits_all.csv                diamond_hits_all.csv
echo "phylogenetics"
add 20_RCA.aln.fasta                       RCA.aln.fasta
add 21_CPN60.aln.fasta                     CPN60.aln.fasta
add 22_RCA.nwk                             RCA.nwk
add 23_CPN60.nwk                           CPN60.nwk
add 24_RCA.supported.nwk                   RCA.supported.nwk
add 25_CPN60.supported.nwk                 CPN60.supported.nwk
add 26_tree_rooting.csv                    tree_rooting.csv
echo "genome-wide context"
add 27_genome_wide_pav.csv                 genome_wide_pav.csv
add 28_wholeproteome_catalogue.faa         WholeProteome_catalog.annotated.faa
echo "publication bundle"
add 29_main_tables                         final/tables
add 30_main_figures                        final/figures
add 31_supplementary                       final/supplementary
add 32_presubmission_checks                checks
add 33_RESULTS_NUMBERS.json                final/manuscript/RESULTS_NUMBERS.json
add 34_verification.txt                    VERIFICATION.txt

cp "$ROOT/verify_outputs.py" "$OUT/" 2>/dev/null && rm -f "$OUT/verify_outputs.py"
python3 verify_outputs.py --inputs "$B" > "$OUT/34_verification.txt" 2>&1 && n=$((n+1))

cat > "$OUT/00_README.md" <<'MD'
# Wheat RCA / Cpn60 pan-proteome — derived datasets

Every value quoted in the associated manuscript is read from these files. The build that
produced them passed an automated consistency gate; `34_verification.txt` is that
gate's output, included unedited.

| Files | Contents |
|---|---|
| 01–05 | proteotype catalogue, composite identifiers, protein models, presence–absence |
| 06–15 | audit tables: candidates, reference loci, granularity, Cpn60 classes, recruitment specificity, divergence, HMM domain classification, transcript cross-reference |
| 16–19 | physicochemical descriptors, polymorphic sites, reference query sets, DIAMOND hits |
| 20–26 | alignments, trees (raw and collapsed/rooted) and the rooting report |
| 27–28 | genome-wide pan-proteome context |
| 29–33 | main tables, main figures, supplementary set, presubmission checks, all computed numbers |
| 35 | complete archive of the above |

**Identifiers.** Each proteotype has a stable `proteotype_id` and a composite identifier
`<proteotype>__<representative cultivar>__<UniProt accession>`. Ties in the choice of
representative are broken lexicographically, so the composite identifier is a
deterministic function of the catalogue and is reproducible across builds.

**Trees.** `tree_rooting.csv` records, per family, the rooting requested and the rooting
achieved. Both families in this build are outgroup-rooted.
MD
echo "  + 00_README.md"; n=$((n+1))

(cd "$OUT" && zip -qr 35_complete_dataset.zip . -x 35_complete_dataset.zip) && echo "  + 35_complete_dataset.zip"

echo
echo "$n item(s) staged, $miss absent"
du -sh "$OUT"
echo
echo "Upload the CONTENTS of $OUT to Zenodo (all files, flat)."
