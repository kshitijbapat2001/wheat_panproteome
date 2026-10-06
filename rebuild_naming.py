#!/usr/bin/env python3
"""Rebuild proteotype_naming.csv with deterministic composite identifiers and report"""
import argparse, os, shutil, sys, glob, csv, datetime

def find(root, name):
    hits = glob.glob(os.path.join(root, "**", name), recursive=True)
    return sorted(hits, key=len)[0] if hits else None

def load(path):
    if not path or not os.path.exists(path): return {}
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["proteotype_id"]: r for r in csv.DictReader(fh)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True, help="build directory, e.g. outputs/latest")
    ap.add_argument("--code", default=".", help="directory holding hpc_pipeline.py")
    ap.add_argument("--dry-run", action="store_true", help="report changes, write nothing")
    a = ap.parse_args()

    need = ["proteotypes.csv", "members.csv", "diamond_hits_all.csv"]
    found = {n: find(a.inputs, n) for n in need}
    for n, p in found.items():
        if not p: sys.exit(f"required file not found under {a.inputs!r}: {n}")
        print(f"[naming] input  : {p}")
    work = os.path.dirname(found["proteotypes.csv"])

    sys.path.insert(0, os.path.abspath(a.code))
    try:
        import hpc_pipeline as HP
    except Exception as e:
        sys.exit(f"cannot import hpc_pipeline from {a.code!r}: {type(e).__name__}: {e}")

    src = open(os.path.join(a.code, "hpc_pipeline.py"), encoding="utf-8").read()
    if 'sort_values(\n                ["evalue_to_ref", "cultivar", "source_id"]' not in src:
        print("\n[naming] WARNING: hpc_pipeline.py does not carry the determinism fix.")
        print("[naming]          run fix_composite_id_stability.py first, or names stay unstable.\n")

    old_path = find(a.inputs, "proteotype_naming.csv")
    old = load(old_path)
    if old_path: print(f"[naming] compare: {old_path}  ({len(old)} rows)")

    out = os.path.join(work, "proteotype_naming.NEW.csv")
    HP.build_proteotype_naming(work, out)
    new = load(out)
    print(f"[naming] wrote  : {out}  ({len(new)} rows)\n")

    if not old:
        print("[naming] no previous table to compare against."); return 0

    changed, gone, added = [], [], []
    for pid, r in new.items():
        if pid not in old: added.append(pid)
        elif old[pid].get("composite_id") != r.get("composite_id"):
            changed.append((pid, old[pid]["composite_id"], r["composite_id"]))
    gone = [p for p in old if p not in new]

    print(f"{'proteotype':<15}{'was':<42}{'now':<42}")
    for pid, o, n in sorted(changed):
        print(f"{pid:<15}{o:<42}{n:<42}")
    print()
    print(f"  changed : {len(changed)} of {len(new)}")
    print(f"  added   : {len(added)}{'  ' + ', '.join(added[:6]) if added else ''}")
    print(f"  removed : {len(gone)}{'  ' + ', '.join(gone[:6]) if gone else ''}")

    if changed:
        rep = os.path.join(work, "composite_id_changes.csv")
        if not a.dry_run:
            with open(rep, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh); w.writerow(["proteotype_id", "old_composite_id", "new_composite_id"])
                w.writerows(sorted(changed))
            print(f"\n  change log -> {rep}")
        print("\n  These artefacts embed composite identifiers and must be regenerated:")
        print("    - Table 5 (composite_id column)")
        print("    - Figure 2 and Supplementary Figures S1/S2 (row labels)")
        print("    - the annotated FASTA catalogue (sequence headers)")
        print("    - Supplementary Table S7 (proteotype_naming.csv)")

    if a.dry_run:
        os.remove(out); print("\n[naming] --dry-run: NEW file removed, nothing replaced.")
        return 0
    if old_path:
        bak = old_path + f".bak_{datetime.datetime.now():%Y%m%d}"
        if not os.path.exists(bak): shutil.copy(old_path, bak); print(f"[naming] backup -> {bak}")
        shutil.copy(out, old_path); print(f"[naming] replaced -> {old_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
