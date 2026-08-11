#!/usr/bin/env python3
"""Apply the v3 figure-size fixes to files already in this directory.

Three changes, all idempotent:
  1. verify_outputs.py  - lift Pillow's decompression-bomb limit, report oversized
                          figures as a normal failure instead of raising
  2. hpc_pipeline.py    - clamp savefig dpi so no figure exceeds 60 Mpx (never < 300 dpi)
                        - cap the content-scaled heatmap height at 30 inches
  3. make_all_figures.py - the same savefig clamp

Run:  python apply_v3_fixes.py
"""
import ast, os, shutil, sys

GUARD = '''
# --- figure size guard ---------------------------------------------------------------
# Figure size scales with content here (a 141-row heatmap is ~58 inches tall). At 600 dpi
# that is 189 megapixels: slow to write, impossible to open, and Pillow refuses to read
# it at all ("decompression bomb"). Clamp the effective dpi so no figure exceeds
# MAX_FIG_PIXELS, never going below 300, which is the usual journal minimum.
MAX_FIG_PIXELS = 60e6
MIN_FIG_DPI = 300
try:
    import matplotlib.figure as _mfig
    _orig_savefig = _mfig.Figure.savefig

    def _clamped_savefig(self, *args, **kw):
        w, h = self.get_size_inches()
        want = kw.get("dpi") or plt.rcParams.get("savefig.dpi", 300)
        if not isinstance(want, (int, float)):
            want = 300.0
        cap = (MAX_FIG_PIXELS / max(float(w) * float(h), 1e-6)) ** 0.5
        kw["dpi"] = max(MIN_FIG_DPI, min(float(want), cap))
        return _orig_savefig(self, *args, **kw)

    _mfig.Figure.savefig = _clamped_savefig
except Exception:  # noqa: BLE001
    pass
'''

OLD_FIG = '''    # 6 -- figure resolution
    try:
        from PIL import Image
        low = []
        for p in glob.glob(os.path.join(L, "**", "*.png"), recursive=True):
            w, _ = Image.open(p).size
            if w / 6.5 < a.min_dpi: low.append((os.path.basename(p), round(w/6.5)))
        (r.ok if not low else r.bad)(
            f"all figures >= {a.min_dpi:.0f} dpi at 6.5in" if not low
            else f"low-resolution figures: {low}")
    except ImportError:
        r.warn("Pillow not installed; figure DPI not checked")'''

NEW_FIG = '''    # 6 -- figure resolution and size
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
        r.warn("Pillow not installed; figure DPI not checked")'''

OLD_HM = '''    fig, ax = plt.subplots(figsize=(max(6, 0.6 * mat.shape[1] + 3),
                                    max(3, 0.4 * mat.shape[0] + 2)))'''
NEW_HM = '''    # Height is capped: beyond ~30 inches the rows are unreadable anyway, and the figure
    # becomes too large to open. Tall matrices belong in the supplementary set.
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * mat.shape[1] + 3),
                                    min(30.0, max(3, 0.4 * mat.shape[0] + 2))))'''

ANCHORS = {
    "hpc_pipeline.py": 'plt.rcParams.update({"figure.dpi": 600, "savefig.dpi": 600})\n',
    "make_all_figures.py": '    "axes.edgecolor": "#333333", "figure.dpi": 300, "savefig.dpi": 300,\n})\n',
}

# shrink_oversized_figures.py, base64 so its own triple quotes cannot break this file
SHRINK_B64 = """
IyEvdXNyL2Jpbi9lbnYgcHl0aG9uMwoiIiJEb3duc2FtcGxlIGZpZ3VyZXMgdGhhdCBhcmUgdG9vIGxhcmdlIHRvIG9wZW4s
IGtlZXBpbmcgdGhlbSBhYm92ZSAzMDAgZHBpLgoKQSBjb250ZW50LXNjYWxlZCBoZWF0bWFwIHJlbmRlcmVkIGF0IDYwMCBk
cGkgY2FuIHJlYWNoIH4xOTAgbWVnYXBpeGVscy4gUGlsbG93IHJlZnVzZXMKdG8gcmVhZCBzdWNoIGEgZmlsZSwgam91cm5h
bHMgd2lsbCBub3QgYWNjZXB0IGl0LCBhbmQgaXQgaXMgdW5yZWFkYWJsZSBhbnl3YXkuIFRoaXMKcmVzY2FsZXMgYW55IFBO
RyBhYm92ZSAtLW1heC1waXhlbHMgc28gdGhhdCBpdCBzdGF5cyBhdCBvciBhYm92ZSAtLW1pbi1kcGkgYWNyb3NzIGEKNi41
LWluY2ggY29sdW1uLCBhbmQga2VlcHMgdGhlIG9yaWdpbmFsIGFsb25nc2lkZSBhcyAub3JpZy5wbmcuCgpPbmx5IG5lZWRl
ZCBmb3IgYnVpbGRzIHByb2R1Y2VkIGJlZm9yZSB0aGUgc2F2ZWZpZyBjbGFtcCB3YXMgYWRkZWQ7IG5ldyBydW5zIGNhbm5v
dApjcmVhdGUgc3VjaCBmaWxlcy4KClVzYWdlCiAgICBweXRob24gc2hyaW5rX292ZXJzaXplZF9maWd1cmVzLnB5IC0taW5w
dXRzIG91dHB1dHMvbGF0ZXN0CiAgICBweXRob24gc2hyaW5rX292ZXJzaXplZF9maWd1cmVzLnB5IC0taW5wdXRzIDxidWls
ZD4gLS1kcnktcnVuCiIiIgppbXBvcnQgYXJncGFyc2UsIGdsb2IsIG9zLCBzeXMKCmRlZiBtYWluKCk6CiAgICBhcCA9IGFy
Z3BhcnNlLkFyZ3VtZW50UGFyc2VyKCkKICAgIGFwLmFkZF9hcmd1bWVudCgiLS1pbnB1dHMiLCByZXF1aXJlZD1UcnVlKQog
ICAgYXAuYWRkX2FyZ3VtZW50KCItLW1heC1waXhlbHMiLCB0eXBlPWZsb2F0LCBkZWZhdWx0PTYwZTYpCiAgICBhcC5hZGRf
YXJndW1lbnQoIi0tbWluLWRwaSIsIHR5cGU9ZmxvYXQsIGRlZmF1bHQ9MzAwKQogICAgYXAuYWRkX2FyZ3VtZW50KCItLXdp
ZHRoLWluIiwgdHlwZT1mbG9hdCwgZGVmYXVsdD02LjUpCiAgICBhcC5hZGRfYXJndW1lbnQoIi0tZHJ5LXJ1biIsIGFjdGlv
bj0ic3RvcmVfdHJ1ZSIpCiAgICBhID0gYXAucGFyc2VfYXJncygpCiAgICB0cnk6CiAgICAgICAgZnJvbSBQSUwgaW1wb3J0
IEltYWdlCiAgICBleGNlcHQgSW1wb3J0RXJyb3I6CiAgICAgICAgc3lzLmV4aXQoIlBpbGxvdyByZXF1aXJlZDogIGNvbmRh
IGluc3RhbGwgcGlsbG93IikKICAgIEltYWdlLk1BWF9JTUFHRV9QSVhFTFMgPSBOb25lCgogICAgcG5ncyA9IHNvcnRlZChn
bG9iLmdsb2Iob3MucGF0aC5qb2luKGEuaW5wdXRzLCAiKioiLCAiKi5wbmciKSwgcmVjdXJzaXZlPVRydWUpKQogICAgcG5n
cyA9IFtwIGZvciBwIGluIHBuZ3MgaWYgbm90IHAuZW5kc3dpdGgoIi5vcmlnLnBuZyIpXQogICAgaWYgbm90IHBuZ3M6CiAg
ICAgICAgcHJpbnQoZiJubyBQTkdzIHVuZGVyIHthLmlucHV0c30iKTsgcmV0dXJuIDAKICAgIHByaW50KGYie2xlbihwbmdz
KX0gZmlndXJlKHMpIHVuZGVyIHthLmlucHV0c31cbiIpCiAgICBuID0gMAogICAgZm9yIHAgaW4gcG5nczoKICAgICAgICB0
cnk6CiAgICAgICAgICAgIGltID0gSW1hZ2Uub3BlbihwKTsgdywgaCA9IGltLnNpemUKICAgICAgICBleGNlcHQgRXhjZXB0
aW9uIGFzIGU6ICAjIG5vcWE6IEJMRTAwMQogICAgICAgICAgICBwcmludChmIiAgVU5SRUFEQUJMRSB7b3MucGF0aC5iYXNl
bmFtZShwKX06IHt0eXBlKGUpLl9fbmFtZV9ffSIpOyBjb250aW51ZQogICAgICAgIGlmIHcgKiBoIDw9IGEubWF4X3BpeGVs
czoKICAgICAgICAgICAgY29udGludWUKICAgICAgICBzY2FsZSA9IChhLm1heF9waXhlbHMgLyAodyAqIGgpKSAqKiAwLjUK
ICAgICAgICBudywgbmggPSBtYXgoMSwgaW50KHcgKiBzY2FsZSkpLCBtYXgoMSwgaW50KGggKiBzY2FsZSkpCiAgICAgICAg
ZWZmID0gbncgLyBhLndpZHRoX2luCiAgICAgICAgbm90ZSA9ICIiIGlmIGVmZiA+PSBhLm1pbl9kcGkgZWxzZSBmIiAgKG9u
bHkge2VmZjouMGZ9IGRwaSBhdCB7YS53aWR0aF9pbn1pbiAtIHZlcnkgdGFsbCBmaWd1cmUpIgogICAgICAgIHByaW50KGYi
ICB7b3MucGF0aC5iYXNlbmFtZShwKX0iKQogICAgICAgIHByaW50KGYiICAgIHt3fXh7aH0gKHt3KmgvMWU2Oi4wZn0gTXB4
KSAtPiB7bnd9eHtuaH0gKHtudypuaC8xZTY6LjBmfSBNcHgpe25vdGV9IikKICAgICAgICBpZiBhLmRyeV9ydW46CiAgICAg
ICAgICAgIG4gKz0gMTsgY29udGludWUKICAgICAgICBvcmlnID0gcFs6LTRdICsgIi5vcmlnLnBuZyIKICAgICAgICBpZiBu
b3Qgb3MucGF0aC5leGlzdHMob3JpZyk6CiAgICAgICAgICAgIG9zLnJlbmFtZShwLCBvcmlnKQogICAgICAgICAgICBpbSA9
IEltYWdlLm9wZW4ob3JpZykKICAgICAgICBpbS5jb252ZXJ0KCJSR0IiKS5yZXNpemUoKG53LCBuaCksIEltYWdlLkxBTkNa
T1MpLnNhdmUocCwgb3B0aW1pemU9VHJ1ZSkKICAgICAgICBuICs9IDEKICAgIHByaW50KCkKICAgIGlmIG4gPT0gMDoKICAg
ICAgICBwcmludCgiICBub3RoaW5nIGFib3ZlIHRoZSBsaW1pdDsgbm8gY2hhbmdlLiIpCiAgICBlbGlmIGEuZHJ5X3J1bjoK
ICAgICAgICBwcmludChmIiAgLS1kcnktcnVuOiB7bn0gZmlndXJlKHMpIHdvdWxkIGJlIHJlc2NhbGVkLiIpCiAgICBlbHNl
OgogICAgICAgIHByaW50KGYiICB7bn0gZmlndXJlKHMpIHJlc2NhbGVkOyBvcmlnaW5hbHMga2VwdCBhcyAqLm9yaWcucG5n
IikKICAgIHJldHVybiAwCgppZiBfX25hbWVfXyA9PSAiX19tYWluX18iOgogICAgc3lzLmV4aXQobWFpbigpKQo=
"""


def write_shrink():
    """Always (re)write shrink_oversized_figures.py so a hand-pasted copy cannot differ."""
    import base64
    want = base64.b64decode("".join(SHRINK_B64.split())).decode("utf-8")
    p = "shrink_oversized_figures.py"
    cur = open(p, encoding="utf-8").read() if os.path.exists(p) else ""
    if cur == want:
        print(f"  ok   {p}: already correct"); return True
    try:
        ast.parse(want)
    except SyntaxError as e:
        print(f"  FAIL {p}: embedded copy is broken ({e})"); return False
    if cur and not os.path.exists(p + ".pre_v3"):
        shutil.copy(p, p + ".pre_v3")
    open(p, "w", encoding="utf-8").write(want)
    print(f"  OK   {p}: written from the canonical copy")
    return True


def edit(fname, fn):
    if not os.path.exists(fname):
        print(f"  MISS {fname} not found here"); return False
    src = open(fname, encoding="utf-8").read()
    new = fn(src)
    if new is None:
        return False
    if new == src:
        print(f"  ok   {fname}: already up to date"); return True
    try:
        ast.parse(new)
    except SyntaxError as e:
        print(f"  FAIL {fname}: patch would break syntax at line {e.lineno}"); return False
    if not os.path.exists(fname + ".pre_v3"):
        shutil.copy(fname, fname + ".pre_v3")
    open(fname, "w", encoding="utf-8").write(new)
    print(f"  OK   {fname}: patched (backup {fname}.pre_v3)")
    return True

def patch_verify(s):
    if "MAX_IMAGE_PIXELS" in s: return s
    if OLD_FIG not in s:
        print("  MISS verify_outputs.py: figure block not recognised"); return None
    return s.replace(OLD_FIG, NEW_FIG, 1)

def patch_guard(name):
    def f(s):
        out = s
        if "MAX_FIG_PIXELS" not in out:
            a = ANCHORS[name]
            if a not in out:
                print(f"  MISS {name}: rcParams anchor not found"); return None
            out = out.replace(a, a + GUARD, 1)
        if name == "hpc_pipeline.py" and OLD_HM in out:
            out = out.replace(OLD_HM, NEW_HM, 1)
        return out
    return f

def main():
    here = os.path.dirname(os.path.abspath(__file__)) or "."
    os.chdir(here)
    print(f"applying v3 figure-size fixes in {here}\n")
    ok = True
    ok &= edit("verify_outputs.py", patch_verify)
    ok &= edit("hpc_pipeline.py", patch_guard("hpc_pipeline.py"))
    ok &= edit("make_all_figures.py", patch_guard("make_all_figures.py"))
    ok &= write_shrink()
    print()
    for f, marker in (("verify_outputs.py", "MAX_IMAGE_PIXELS"),
                      ("hpc_pipeline.py", "MAX_FIG_PIXELS"),
                      ("make_all_figures.py", "MAX_FIG_PIXELS")):
        present = os.path.exists(f) and marker in open(f, encoding="utf-8").read()
        print(f"  {'OK  ' if present else 'FAIL'} {f} carries {marker}")
        ok &= present
    cap = os.path.exists("hpc_pipeline.py") and "min(30.0, max(3, 0.4" in open("hpc_pipeline.py", encoding="utf-8").read()
    print(f"  {'OK  ' if cap else 'FAIL'} hpc_pipeline.py heatmap height capped")
    print("\n" + ("all fixes in place." if ok and cap else "SOMETHING IS MISSING - see above."))
    return 0 if (ok and cap) else 1

if __name__ == "__main__":
    sys.exit(main())
