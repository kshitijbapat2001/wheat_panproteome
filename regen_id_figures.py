#!/usr/bin/env python3
"""Regenerate only the figures whose row labels embed a composite identifier."""
import argparse, os, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--suppdir", required=True)
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--cpn60-rows", type=int, default=30,
                    help="rows shown in Figure 2b; must match the caption (default 30)")
    ap.add_argument("--code", default=".")
    ap.add_argument("--label", choices=["composite", "id"], default="composite",
                    help="composite reproduces the submitted figures")
    a = ap.parse_args()

    sys.path.insert(0, os.path.abspath(a.code))
    import matplotlib; matplotlib.use("Agg")
    try:
        import make_all_figures as M
    except ImportError as e:
        sys.exit(f"cannot import make_all_figures from {a.code!r}: {e}")
    import pandas as pd

    src = open(os.path.join(a.code, "make_all_figures.py"), encoding="utf-8").read()
    if 'sort_values(\n                ["evalue_to_ref", "cultivar", "source_id"]' not in src:
        print("[figs] WARNING: make_all_figures.py lacks the determinism fix; labels may")
        print("[figs]          not match the corrected proteotype_naming.csv.\n")

    M.plt.rcParams.update({"figure.dpi": a.dpi, "savefig.dpi": a.dpi})
    os.makedirs(a.outdir, exist_ok=True); os.makedirs(a.suppdir, exist_ok=True)

    def need(name):
        p = M._find_input(a.inputs, name)
        if not p: sys.exit(f"required input not found under {a.inputs!r}: {name}")
        return p

    pt = pd.read_csv(need("proteotypes.csv"))
    pt["novel_vs_reference"] = pt["novel_vs_reference"].astype(str).str.lower().isin(["true", "1"])
    pav = pd.read_csv(need("pav_matrix.csv"))
    cults = M.cultivar_columns(pav)
    print(f"[figs] {len(pt)} proteotypes x {len(cults)} cultivars, {a.dpi} dpi")

    labels = None
    naming_dir = os.path.dirname(need("members.csv"))
    nm = M.build_proteotype_naming(naming_dir, os.path.join(a.outdir, "proteotype_naming.csv"))
    if a.label == "composite":
        labels = M.composite_label_map(nm)
        print(f"[figs] composite labels built for {len(labels)} proteotypes")

    from PIL import Image
    jobs = []
    f2 = os.path.join(a.outdir, "Figure2_novel_pav.png")
    M.plot_novel_pav(pt, pav, cults, f2, cpn60_rows=a.cpn60_rows, labels=labels)
    jobs.append(("Figure2_novel_pav.png", "manuscript Figure 2", f2))
    for fam, fname, role in [("RCA", "FigureS1_RCA_PAV_heatmap.png", "Supplementary Figure S1"),
                             ("CPN60", "FigureS2_CPN60_PAV_heatmap.png", "Supplementary Figure S2")]:
        p = os.path.join(a.suppdir, fname)
        M.plot_pav_heatmap(pt, pav, cults, fam, p, labels=labels)
        jobs.append((fname, role, p))

    print()
    print(f"{'file':<34}{'as':<26}{'pixels':>14}{'dpi@6.5in':>11}   note")
    for fname, role, p in jobs:
        w, h = Image.open(p).size
        eff = w / 6.5
        note = ""
        if h / max(w, 1) > 2.9:
            note = "very tall - place landscape or split"
        print(f"{fname:<34}{role:<26}{w:>6}x{h:<7}{eff:>10.0f}   {note}")
    print(f"\n[figs] wrote to {os.path.abspath(a.outdir)} and {os.path.abspath(a.suppdir)}")
    print("[figs] Figure 2b shows the %d most prevalent Cpn60 rows - the caption says 30." % a.cpn60_rows)
    return 0

if __name__ == "__main__":
    sys.exit(main())
