#!/usr/bin/env python3
"""Re-render the three low-resolution genome-wide figures at publication DPI.

Figure 7 (pan-genome overview), Supplementary Figure S3 (prevalence spectrum) and
Supplementary Figure S4 (per-cultivar composition) in the submitted manuscript were
taken from hpc_pipeline.py, whose copies of these plotting functions set no DPI and
therefore save at matplotlib's default 100. make_all_figures.py already contains
correct versions at 300 DPI with (a)/(b) panel letters; this script calls those,
at a DPI you choose, without re-running the pipeline.

Usage
    python regen_lowres_figures.py --inputs outputs/latest --outdir outputs/latest/final/figures
    python regen_lowres_figures.py --inputs outputs/latest --outdir /tmp/figs --dpi 600
"""
import argparse, os, sys, glob

def find(root, name):
    hits = glob.glob(os.path.join(root, "**", name), recursive=True)
    return sorted(hits, key=len)[0] if hits else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True, help="build directory, e.g. outputs/latest")
    ap.add_argument("--outdir", required=True, help="where to write the PNGs")
    ap.add_argument("--dpi", type=int, default=600,
                    help="600 for line art (default); 300 is the usual minimum")
    ap.add_argument("--code", default=".", help="directory holding make_all_figures.py")
    a = ap.parse_args()

    sys.path.insert(0, os.path.abspath(a.code))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import make_all_figures as M
    except ImportError as e:
        sys.exit(f"cannot import make_all_figures from {a.code!r}: {e}")
    import pandas as pd

    # Raise DPI globally. make_all_figures sets 300 at import; override after it.
    M.plt.rcParams.update({"figure.dpi": a.dpi, "savefig.dpi": a.dpi})

    gw_path = find(a.inputs, "genome_wide_pav.csv")
    if not gw_path:
        sys.exit(f"genome_wide_pav.csv not found under {a.inputs!r}. It is written only when "
                 "whole_proteome.enabled is true in config.yaml.")
    print(f"[regen] table : {gw_path}")
    gw = pd.read_csv(gw_path)
    cults = [c for c in gw.columns if c not in ("cluster_id", "family", "novel_vs_reference")]
    print(f"[regen] {len(gw)} clusters x {len(cults)} cultivars, target {a.dpi} dpi")

    os.makedirs(a.outdir, exist_ok=True)
    jobs = [
        ("Figure5_pangenome_overview.png",      "manuscript Figure 7",
         lambda p: M.plot_pangenome_overview(gw, cults, p)),
        ("FigureS5_prevalence_spectrum.png",    "Supplementary Figure S3",
         lambda p: M.plot_prevalence_hist(gw, cults, p)),
        ("FigureS6_per_cultivar_composition.png", "Supplementary Figure S4",
         lambda p: M.plot_genomewide_per_cultivar(gw, cults, p)),
    ]
    from PIL import Image
    print()
    print(f"{'file':<42}{'as':<26}{'pixels':>14}{'dpi @6.5in':>12}")
    ok = True
    for fname, role, fn in jobs:
        out = os.path.join(a.outdir, fname)
        fn(out)
        w, h = Image.open(out).size
        eff = w / 6.5
        flag = "OK" if eff >= 300 else "STILL LOW"
        if eff < 300: ok = False
        print(f"{fname:<42}{role:<26}{w:>6}x{h:<7}{eff:>9.0f}  {flag}")
    print()
    print(f"[regen] wrote to {os.path.abspath(a.outdir)}")
    if ok:
        print("[regen] all three now exceed 300 dpi at 6.5-inch print width.")
        print("[regen] Figure 7 also carries (a)/(b) panel letters, which the caption requires.")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
