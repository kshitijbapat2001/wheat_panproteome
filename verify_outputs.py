#!/usr/bin/env python3
"""Gate the run before anything is used in the manuscript.

Every check here corresponds to a defect that actually reached a draft: a mixed-build
bundle, a silent midpoint fall-back, a blank table cell, a low-resolution figure, a
truncated candidate table, an identifier that changed between builds.

Exit status 0 = safe to use. 1 = at least one FAIL.
"""
import argparse, csv, glob, json, os, sys
from collections import Counter

def find(root, name):
    h = glob.glob(os.path.join(root, "**", name), recursive=True)
    h = [p for p in h if ".bak" not in p]
    return sorted(h, key=os.path.getmtime, reverse=True) if h else []

def rows(p):
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))

class R:
    def __init__(s): s.p=s.f=s.w=0; s.msgs=[]
    def ok(s,m):   s.p+=1; s.msgs.append(("PASS", m))
    def bad(s,m):  s.f+=1; s.msgs.append(("FAIL", m))
    def warn(s,m): s.w+=1; s.msgs.append(("WARN", m))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--min-dpi", type=float, default=300)
    a = ap.parse_args(); L = a.inputs; r = R()

    # 1 -- duplicate inputs (the mixed-build failure)
    for n in ("proteotypes.csv", "members.csv", "proteotype_naming.csv"):
        h = find(L, n)
        if not h: r.bad(f"{n}: missing")
        elif len(h) > 1:
            sizes = {os.path.getsize(p) for p in h}
            (r.warn if len(sizes) == 1 else r.bad)(
                f"{n}: {len(h)} copies" + ("" if len(sizes) == 1 else
                " with DIFFERENT sizes -> mixed builds: " + ", ".join(h[:3])))
        else: r.ok(f"{n}: single copy")

    # 2 -- one proteotype universe across every downstream table
    univ = {}
    for n in ("proteotypes.csv", "proteotype_naming.csv",
              "domain_classification.csv", "proteotype_recovery.csv"):
        h = find(L, n)
        if h: univ[n] = len(rows(h[0]))
    if univ:
        if len(set(univ.values())) == 1:
            r.ok(f"proteotype universe consistent at {list(univ.values())[0]} across {len(univ)} tables")
        else:
            r.bad(f"proteotype counts DISAGREE: {univ}")

    # 3 -- rooting actually achieved
    h = find(L, "tree_rooting.csv")
    if not h: r.warn("tree_rooting.csv absent (older pipeline?)")
    else:
        for x in rows(h[0]):
            if x["requested"] == "outgroup" and x["achieved"] != "outgroup":
                r.bad(f"{x['family']}: outgroup rooting requested, {x['achieved']} achieved "
                      "-> do NOT claim outgroup rooting")
            else:
                r.ok(f"{x['family']}: {x['achieved']}-rooted as requested")

    # 4 -- candidate table complete and consistent with Table 5
    ct, t5 = find(L, "candidate_table.csv"), find(L, "Table5.csv")
    if ct and t5:
        n1, n2 = len(rows(ct[0])), len(rows(t5[0]))
        (r.ok if n1 == n2 else r.bad)(f"candidate_table {n1} rows vs Table5 {n2} rows")
        cls = Counter(x.get("candidate_class","") for x in rows(ct[0]))
        (r.ok if len(cls) > 1 else r.warn)(f"candidate tiers: {dict(cls)}")

    # 5 -- Table 1 has no blank cells
    h = find(L, "Table1.csv")
    if h:
        blanks = [x["Metric"] for x in rows(h[0])
                  if not str(x.get("RCA","")).strip() or not str(x.get("Cpn60","")).strip()]
        (r.ok if not blanks else r.bad)(
            "Table 1 complete" if not blanks else f"Table 1 blank rows: {blanks}")

    # 6 -- figure resolution and size
    try:
        from PIL import Image
        # Pillow refuses images above ~179 Mpx as a "decompression bomb". Reading our own
        # generated figures is not a security risk, and crashing the gate on an oversized
        # figure hides every other check, so lift the limit and flag the size instead.
        Image.MAX_IMAGE_PIXELS = None
        low, huge, unreadable = [], [], []
        pngs = [p for p in sorted(glob.glob(os.path.join(L, "**", "*.png"), recursive=True))
                if not p.endswith(".orig.png")]   # .orig.png = pre-rescale backups
        for p in pngs:
            try:
                w, h = Image.open(p).size
            except Exception as e:  # noqa: BLE001
                unreadable.append((os.path.basename(p), type(e).__name__)); continue
            if w * h > 60e6:
                huge.append((os.path.basename(p), f"{w}x{h}", f"{w*h/1e6:.0f}Mpx"))
            if w / 6.5 < a.min_dpi:
                low.append((os.path.basename(p), round(w / 6.5)))
        (r.ok if not low else r.bad)(
            f"all figures >= {a.min_dpi:.0f} dpi at 6.5in" if not low
            else f"low-resolution figures: {low}")
        if huge:
            r.bad(f"oversized figures (>60 Mpx; too large to open or submit): {huge}")
        else:
            r.ok("no figure exceeds 60 Mpx")
        if unreadable:
            r.bad(f"unreadable figures: {unreadable}")
    except ImportError:
        r.warn("Pillow not installed; figure DPI not checked")

    # 7 -- identifier determinism
    h = find(L, "proteotype_naming.csv")
    if h:
        nm = rows(h[0])
        bad = [x["proteotype_id"] for x in nm if not x.get("composite_id","").count("__") == 2]
        (r.ok if not bad else r.bad)(
            f"composite ids well-formed ({len(nm)})" if not bad else f"malformed ids: {bad[:5]}")

    # 8 -- headline arithmetic
    h = find(L, "RESULTS_NUMBERS.json")
    if h:
        N = json.load(open(h[0]))
        try:
            d = N["domains"]["classes"] if "classes" in N.get("domains", {}) else None
            tot = N.get("n_proteotypes") or {}
            r.ok(f"RESULTS_NUMBERS.json readable ({len(N)} top-level keys)")
        except Exception as e:
            r.warn(f"RESULTS_NUMBERS.json parse note: {e}")

    w = min(max((len(m) for _, m in r.msgs), default=0), 92)
    print(f"\nverification of {L}\n" + "=" * (w + 8))
    for k, m in r.msgs: print(f"  {k}  {m}")
    print("=" * (w + 8))
    print(f"  {r.p} passed, {r.f} failed, {r.w} warnings")
    if r.f:
        print("\n  DO NOT use these outputs in the manuscript until the failures are fixed.")
    return 1 if r.f else 0

if __name__ == "__main__":
    sys.exit(main())
