#!/usr/bin/env bash
# Print the canonical path, date and key numbers for every file needed to edit the
# manuscript, so a stale copy is obvious BEFORE it is uploaded.
#
#   bash whats_fresh.sh                       # newest build under the configured outputs
#   bash whats_fresh.sh outputs/pan_proteome_build_YYYYMMDD_HHMMSS
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
B="${1:-}"
[[ -z "$B" ]] && B="$(ls -1dt outputs/pan_proteome_build_* 2>/dev/null | head -1)"
[[ -n "$B" && -d "$B" ]] || { echo "no build found; pass one explicitly"; exit 1; }
B="$(readlink -f "$B")"
echo "build: $B"
echo

show () {  # show <label> <relative path> [awk-summary]
  local lbl="$1" rel="$2"; shift 2
  local p="$B/$rel"
  if [[ -f "$p" ]]; then
    printf "  %-26s %s  %s\n" "$lbl" "$(date -r "$p" +%Y-%m-%d_%H:%M)" "$rel"
    [[ $# -gt 0 ]] && { printf "    %s\n" "$("$@" "$p")"; }
  else
    printf "  %-26s MISSING            %s\n" "$lbl" "$rel"
  fi
}

conf () { python3 -c "
import csv,collections,sys
d=list(csv.DictReader(open(sys.argv[1])))
c=collections.Counter(r.get('confidence','?') for r in d)
print(f'{len(d)} rows | ' + ', '.join(f'{k}={v}' for k,v in c.most_common()))" "$1"; }
gen ()  { python3 -c "
import json,sys; N=json.load(open(sys.argv[1]))
d=N.get('domains',{})
print(f\"generated={N.get('generated')} | domains n={d.get('n')} conf={d.get('confidence')}\")" "$1"; }
rows () { python3 -c "
import csv,sys; d=list(csv.DictReader(open(sys.argv[1]))); print(f'{len(d)} rows')" "$1"; }
root () { python3 -c "
import csv,sys
for r in csv.DictReader(open(sys.argv[1])):
    print(f\"{r['family']}: requested={r['requested']} achieved={r['achieved']} outgroups=[{r['outgroup_tips']}]\")" "$1"; }

echo "FILES TO SEND"
show "RESULTS_NUMBERS.json"   final/manuscript/RESULTS_NUMBERS.json gen
show "domain_classification"  domains/domain_classification.csv     conf
show "tree_rooting"           tree_rooting.csv                      root
show "candidate_table"        audit/candidate_table.csv             rows
show "proteotypes"            proteotypes.csv                       rows
show "proteotype_naming"      proteotype_naming.csv                 rows
show "members"                members.csv                           rows
show "prevalence_discordance" transcripts/prevalence_discordance.csv rows
show "Table1"                 final/tables/Table1.csv               rows
show "Table5"                 final/tables/Table5.csv               rows
echo
echo "Today is $(date +%Y-%m-%d). Anything dated earlier than the build is stale."
echo "Bundle everything with:  bash run_all.sh collect"
