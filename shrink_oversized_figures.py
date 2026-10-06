#!/usr/bin/env python3
"""Downsample figures that are too large to open, keeping them above 300 dpi."""
import argparse, glob, os, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--max-pixels", type=float, default=60e6)
    ap.add_argument("--min-dpi", type=float, default=300)
    ap.add_argument("--width-in", type=float, default=6.5)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    try:
        from PIL import Image
    except ImportError:
        sys.exit("Pillow required:  conda install pillow")
    Image.MAX_IMAGE_PIXELS = None

    pngs = sorted(glob.glob(os.path.join(a.inputs, "**", "*.png"), recursive=True))
    pngs = [p for p in pngs if not p.endswith(".orig.png")]
    if not pngs:
        print(f"no PNGs under {a.inputs}"); return 0
    print(f"{len(pngs)} figure(s) under {a.inputs}\n")
    n = 0
    for p in pngs:
        try:
            im = Image.open(p); w, h = im.size
        except Exception as e:
            print(f"  UNREADABLE {os.path.basename(p)}: {type(e).__name__}"); continue
        if w * h <= a.max_pixels:
            continue
        scale = (a.max_pixels / (w * h)) ** 0.5
        nw, nh = max(1, int(w * scale)), max(1, int(h * scale))
        eff = nw / a.width_in
        note = "" if eff >= a.min_dpi else f"  (only {eff:.0f} dpi at {a.width_in}in - very tall figure)"
        print(f"  {os.path.basename(p)}")
        print(f"    {w}x{h} ({w*h/1e6:.0f} Mpx) -> {nw}x{nh} ({nw*nh/1e6:.0f} Mpx){note}")
        if a.dry_run:
            n += 1; continue
        orig = p[:-4] + ".orig.png"
        if not os.path.exists(orig):
            os.rename(p, orig)
            im = Image.open(orig)
        im.convert("RGB").resize((nw, nh), Image.LANCZOS).save(p, optimize=True)
        n += 1
    print()
    if n == 0:
        print("  nothing above the limit; no change.")
    elif a.dry_run:
        print(f"  --dry-run: {n} figure(s) would be rescaled.")
    else:
        print(f"  {n} figure(s) rescaled; originals kept as *.orig.png")
    return 0

if __name__ == "__main__":
    sys.exit(main())
