#!/usr/bin/env python3
"""Show exactly which file build_tables.py will pick for each input it needs."""
import os, sys, glob, yaml
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_tables as BT

root = sys.argv[1] if len(sys.argv) > 1 else "outputs/latest"
cfg = yaml.safe_load(open("config.yaml")) or {}
_inp = (cfg.get("project") or {}).get("inputs_dir", "inputs")
BT.EXTRA_SEARCH = [d for d in (_inp, os.path.join(_inp, "references"),
                               os.path.join(root, "..", "..", "inputs"),
                               os.path.join(root, "..", "..", "inputs", "references"))
                   if d and os.path.isdir(d)]
BT.EXCLUDE_DIRS = [os.path.abspath(os.path.join(root, "final"))]
print(f"build        : {os.path.realpath(root)}")
print(f"EXTRA_SEARCH : {BT.EXTRA_SEARCH}\n")
for name in ("reference_query_sets.csv", "proteotype_recovery.csv",
             "domain_classification.csv", "candidate_table.csv",
             "cs_locus_table.csv", "granularity_partition.csv",
             "prevalence_discordance.csv", "physicochemical_features.csv"):
    p = BT._find(root, name)
    print(f"  {name:<32} {'-> ' + p if p else 'NOT FOUND'}")
    if not p:
        hits = glob.glob(os.path.join(root, "**", name), recursive=True)
        hits += glob.glob(os.path.join(root, "..", "..", "**", name), recursive=True)
        seen = {os.path.realpath(h) for h in hits}
        for h in sorted(seen)[:4]:
            print(f"      on disk at: {h}")
        if not seen:
            print("      not anywhere under the build or the project root")
