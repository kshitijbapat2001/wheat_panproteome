#!/usr/bin/env python3
"""Audit figure resolution, either as loose PNGs or as actually placed in a .docx.

For a .docx the effective DPI is what matters: pixel width divided by the width the
image is placed at on the page. A 950-px image is fine as a thumbnail and far too
coarse across a 6.5-inch text column.

Usage
    python check_figure_dpi.py final/figures            # folder of PNGs
    python check_figure_dpi.py Manuscript.docx          # as placed
    python check_figure_dpi.py final/figures --width 3.25 --min-dpi 300
"""
import argparse, io, os, sys, zipfile, glob

EMU = 914400.0
A = '{http://schemas.openxmlformats.org/drawingml/2006/main}'
R = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
WP = '{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}'

def audit_docx(path, min_dpi):
    import docx
    from PIL import Image
    d = docx.Document(path); z = zipfile.ZipFile(path); rels = d.part.rels
    rows = []
    for p in d.paragraphs:
        for blip in p._p.findall('.//' + A + 'blip'):
            tgt = rels[blip.get(R + 'embed')].target_ref
            im = Image.open(io.BytesIO(z.read('word/' + tgt.replace('../', ''))))
            ext = p._p.find('.//' + WP + 'extent')
            if ext is None:
                ext = p._p.find('.//' + A + 'ext')
            w_in = int(ext.get('cx')) / EMU; h_in = int(ext.get('cy')) / EMU
            rows.append((tgt.split('/')[-1], im.size, (w_in, h_in), im.size[0] / w_in))
    return rows

def audit_folder(path, width, min_dpi):
    from PIL import Image
    rows = []
    for f in sorted(glob.glob(os.path.join(path, "*.png")) + glob.glob(os.path.join(path, "*.tif*"))):
        im = Image.open(f)
        rows.append((os.path.basename(f), im.size, (width, width * im.size[1] / im.size[0]),
                     im.size[0] / width))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--width", type=float, default=6.5, help="assumed print width, inches (folder mode)")
    ap.add_argument("--min-dpi", type=float, default=300)
    a = ap.parse_args()

    if a.path.lower().endswith(".docx"):
        rows = audit_docx(a.path, a.min_dpi); mode = "as placed in document"
    elif os.path.isdir(a.path):
        rows = audit_folder(a.path, a.width, a.min_dpi); mode = f"assuming {a.width}in print width"
    else:
        sys.exit("path must be a .docx or a directory")

    if not rows:
        sys.exit("no images found")
    print(f"{a.path}  ({mode}); threshold {a.min_dpi:.0f} dpi\n")
    print(f"{'image':<20}{'pixels':>15}{'placed in':>16}{'eff. dpi':>10}   verdict")
    bad = []
    for name, (pw, ph), (wi, hi), dpi in rows:
        v = "OK" if dpi >= a.min_dpi else ("LOW" if dpi >= a.min_dpi / 2 else "TOO LOW")
        if dpi < a.min_dpi: bad.append((name, dpi))
        # flag extreme aspect ratios: unreadable when squeezed to fit a page
        ar = ""
        if hi > 0 and wi / hi < 0.35: ar = "  <- very tall/narrow, will be illegible"
        print(f"{name:<20}{pw:>7}x{ph:<7}{wi:>6.2f} x {hi:<6.2f}{dpi:>9.0f}   {v}{ar}")
    print()
    if bad:
        print(f"{len(bad)} figure(s) below {a.min_dpi:.0f} dpi:")
        for n, d in bad:
            print(f"  {n}: {d:.0f} dpi -> needs {a.min_dpi/d:.1f}x more pixels")
        return 1
    print("all figures meet the threshold.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
