#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Residue-level alignment views of every proteotype, for protein engineering."""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

CLASSES = [("hydrophobic (AILMV)", "AILMV", "#2a78d6"),
           ("glycine (G)", "G", "#d95926"),
           ("polar (STNQ)", "STNQ", "#199e70"),
           ("proline (P)", "P", "#c98500"),
           ("acidic (DE)", "DE", "#d55181"),
           ("aromatic (FWY)", "FWY", "#008300"),
           ("cysteine (C)", "C", "#4a3aa7"),
           ("basic (KRH)", "KRH", "#e34948")]
COLOUR = {aa: c for _, aas, c in CLASSES for aa in aas}
INK, MUTED, FAINT = "#0b0b0b", "#52514e", "#a9a8a2"

PAGE_W, PAGE_H = 11.69, 8.27
COLS = 100
LINES = 116
FONT = 5.0
MONO = "Courier New"


def _find(root, name):
    for c in (os.path.join(root, name), os.path.join(root, "audit", name),
              os.path.join(root, "domains", name),
              os.path.join(root, "final", "supplementary", "data", name)):
        if os.path.exists(c):
            return c
    hits = glob.glob(os.path.join(root, "**", name), recursive=True)
    return sorted(hits, key=os.path.getmtime, reverse=True)[0] if hits else None


def read_aln(path):
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(buf).upper()))
                name, buf = line[1:].split()[0].split("__")[0], []
            elif line:
                buf.append(line)
    if name:
        out.append((name, "".join(buf).upper()))
    return out


def build_panel(aln, row_ids, ref_id, labels):
    """Subset the alignment to row_ids, drop all-gap columns, number on ref_id."""
    seqs = dict(aln)
    rows = [ref_id] + [r for r in row_ids if r != ref_id]
    keep = [j for j in range(len(seqs[ref_id]))
            if any(seqs[r][j] != "-" for r in rows)]
    cols = {r: "".join(seqs[r][j] for j in keep) for r in rows}
    num, k = [], 0
    for ch in cols[ref_id]:
        if ch != "-":
            k += 1
            num.append(k)
        else:
            num.append(None)
    return {"rows": rows, "seq": cols, "num": num, "labels": labels, "ref": ref_id,
            "features": []}


def col_of(panel, row_id, pos1):
    """Column of residue pos1 (1-based) of row_id within the panel, or None."""
    k = 0
    for j, ch in enumerate(panel["seq"][row_id]):
        if ch != "-":
            k += 1
            if k == pos1:
                return j
    return None


def render(panels, title, pdf_path, png_path, footnote):
    """Lay the panels out block by block and page by page."""
    blocks = []
    for P in panels:
        n = len(P["num"])
        for c0 in range(0, n, COLS):
            need = 2 + sum(2 if t == "domains / regions" else 1
                           for t, _ in P["features_tracks"]) + len(P["rows"]) + 2
            blocks.append((P, c0, need))
    pages, cur, used = [], [], 3
    for b in blocks:
        if used + b[2] > LINES - 3 and cur:
            pages.append(cur); cur, used = [], 3
        cur.append(b); used += b[2]
    if cur:
        pages.append(cur)

    os.makedirs(os.path.dirname(pdf_path), exist_ok=True)
    with PdfPages(pdf_path) as pdf:
        for pi, page in enumerate(pages, 1):
            fig = plt.figure(figsize=(PAGE_W, PAGE_H))
            ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
            ax.set_xlim(-34, COLS + 1)
            ax.set_ylim(LINES, 0)
            ax.axis("off")
            ax.text(-34, 1, f"{title}  (page {pi} of {len(pages)})", fontsize=7,
                    color=INK, weight="bold", va="center", ha="left")
            y = 3
            for P, c0, need in page:
                c1 = min(c0 + COLS, len(P["num"]))
                if c0 == 0:
                    ax.text(-34, y, P["name"], fontsize=5.6, color=INK, weight="bold",
                            va="center", ha="left")
                    y += 1
                else:
                    y += 1
                for tname, feats in P["features_tracks"]:
                    if tname == "domains / regions":
                        y += 1
                    ax.text(-1.2, y, tname, fontsize=FONT - 0.4, color=MUTED, va="center",
                            ha="right", style="italic")
                    for f in feats:
                        a, b = max(f["c0"], c0), min(f["c1"], c1 - 1)
                        if a > b:
                            continue
                        if f.get("kind") == "tick":
                            ax.plot([a - c0], [y], marker="|", ms=3.2, mew=0.6,
                                    color=f.get("colour", INK), ls="none")
                        elif f.get("kind") == "site":
                            ax.plot([a - c0], [y + 0.15], marker="v", ms=2.6,
                                    color=f.get("colour", INK), ls="none")
                            ax.text(a - c0 + 0.8, y - 0.1, f["label"], fontsize=FONT - 0.6,
                                    color=f.get("colour", INK), va="center", ha="left")
                        else:
                            ax.plot([a - c0 - 0.4, b - c0 + 0.4], [y, y], lw=1.6,
                                    color=f.get("colour", INK), solid_capstyle="butt")
                            if f["c0"] >= c0 or a == c0:
                                ax.text(a - c0 - 0.4, y - 0.95, f["label"], fontsize=FONT - 0.3,
                                        color=f.get("colour", INK), va="center", ha="left")
                    y += 1
                for j in range(c0, c1):
                    n = P["num"][j]
                    if n is not None and (n % 10 == 0 or n == 1):
                        ax.text(j - c0, y, str(n), fontsize=FONT - 0.6, color=MUTED,
                                ha="center", va="center")
                ax.text(-1.2, y, "reference position", fontsize=FONT - 0.4, color=MUTED,
                        va="center", ha="right", style="italic")
                y += 1
                ref = P["seq"][P["ref"]]
                for r in P["rows"]:
                    s = P["seq"][r]
                    is_ref = r == P["ref"]
                    ax.text(-1.2, y, P["labels"].get(r, r), fontsize=FONT,
                            color=INK, va="center", ha="right", family=MONO,
                            weight="bold" if is_ref else "normal")
                    dots_x, gap_x = [], []
                    for j in range(c0, c1):
                        ch = s[j]
                        x = j - c0
                        if ch == "-":
                            gap_x.append(x)
                        elif is_ref:
                            ax.text(x, y, ch, fontsize=FONT, color=INK, ha="center",
                                    va="center", family=MONO, weight="bold")
                        elif ch == ref[j]:
                            dots_x.append(x)
                        else:
                            ax.text(x, y, ch, fontsize=FONT, color=COLOUR.get(ch, INK),
                                    ha="center", va="center", family=MONO,
                                    weight="bold")
                    if dots_x:
                        ax.scatter(dots_x, [y] * len(dots_x), s=0.35, c=FAINT, lw=0)
                    if gap_x:
                        ax.scatter(gap_x, [y] * len(gap_x), s=1.4, c=FAINT, marker="_",
                                   lw=0.5)
                    y += 1
                y += 1
            yl = LINES - 1.6
            ax.text(-34, yl - 1.4, "· identical to the reference (top row, bold)    "
                    "– gap    * reference-absent proteoform    (n) cultivars carrying it"
                    "    residue classes:", fontsize=FONT + 0.4, color=MUTED, va="center",
                    ha="left")
            x, step = -34, (COLS + 35) / len(CLASSES)
            for lab, aas, c in CLASSES:
                ax.text(x, yl, aas[0], fontsize=FONT + 1, color=c, weight="bold",
                        family=MONO, va="center", ha="left")
                ax.text(x + 1.6 * COLS / 100, yl, lab, fontsize=FONT + 0.2, color=MUTED,
                        va="center", ha="left")
                x += step
            ax.text(-34, LINES - 0.5, footnote, fontsize=FONT + 0.2, color=MUTED,
                    va="center", ha="left")
            pdf.savefig(fig)
            if pi == 1 and png_path:
                fig.savefig(png_path, dpi=300)
            plt.close(fig)
    return len(pages)


def rca_panel(aln, iso, ref_json):
    ids = [r for r, _ in aln if r.startswith("RCA_PT")]
    order = {"RCA2α": 0, "RCA2β": 1, "RCA1β": 2, "RCA-like": 3}
    iso = iso.set_index("proteotype_id")

    def key(p):
        r = iso.loc[p]
        full = str(r["model_status"]) == "full-length"
        return (order.get(r["isoform"], 9), str(r["established_isoform"]), not full,
                -int(r["n_cultivars"]), p)
    rows = sorted(ids, key=key)
    labels = {}
    for p in rows:
        r = iso.loc[p]
        star = "*" if str(r["novel_vs_reference"]).lower() in ("true", "1") else " "
        var = "" if str(r["model_status"]) == "full-length" else " v"
        labels[p] = f"{p.replace('RCA_', '')} {r['established_isoform']:<14}({int(r['n_cultivars']):>2}){star}{var}"
    alpha, beta = ref_json["alpha_reference"], ref_json["beta_reference"]
    P = build_panel(aln, rows, alpha, labels)
    P["name"] = ("RCA: all 46 proteoforms against the canonical RCA2α precursor "
                 f"({alpha}; residues 1–427 are identical to RCA2β)")
    canon = ref_json.get("canonical", {})
    b1, b2 = canon.get("RCA1β"), canon.get("RCA2β")

    ref_seq = P["seq"][alpha].replace("-", "")
    feats_dom, feats_sites, feats_diff = [], [], []

    def span(label, s1, e1, row=alpha, colour=INK):
        a, b = col_of(P, row, s1), col_of(P, row, e1)
        if a is not None and b is not None:
            feats_dom.append({"label": label, "c0": a, "c1": b, "colour": colour})

    span("transit peptide (1–46)", 1, 46, colour=MUTED)
    m = re.search(r"G[GA]KGQGK[ST]", ref_seq)
    if m:
        span("Walker A", m.start() + 1, m.end(), colour="#2a78d6")
        mb = re.search(r"[ILVMFC]{4}[ND]D", ref_seq[m.end():m.end() + 120])
        if mb:
            span("Walker B", m.end() + mb.start() + 1, m.end() + mb.end(), colour="#2a78d6")
    span("region carrying 8 of 11 thermostabilising substitutions (310–427)",
         310, 427, row=beta, colour="#199e70")
    cte = ref_json["cte_start_1based"]
    span("C-terminal extension (RCA2α only)", cte, len(ref_seq), colour="#4a3aa7")

    sw = col_of(P, beta, ref_json["switch_site_1based_on_beta_reference"])
    if sw is not None:
        feats_sites.append({"label": "M/I switch", "c0": sw,
                            "c1": sw, "kind": "site", "colour": "#e34948"})
    k = col_of(P, alpha, cte)
    if k is not None:
        feats_sites.append({"label": "K428", "c0": k, "c1": k,
                            "kind": "site", "colour": "#e34948"})
    for c in ref_json.get("cte_cysteines_1based", []):
        j = col_of(P, alpha, c)
        if j is not None:
            feats_sites.append({"label": f"C{c}", "c0": j, "c1": j, "kind": "site",
                                "colour": "#4a3aa7"})
    if b1 and b2:
        s1, s2 = P["seq"][b1], P["seq"][b2]
        for j, (x, y) in enumerate(zip(s1, s2)):
            if x != "-" and y != "-" and x != y:
                feats_diff.append({"c0": j, "c1": j, "kind": "tick", "colour": "#d95926"})
    P["features_tracks"] = [("domains / regions", feats_dom),
                            ("functional sites", feats_sites),
                            ("RCA1β ≠ RCA2β", feats_diff)]
    return P


def cpn60_panels(aln, classes, pt):
    ids = [r for r, _ in aln if r.startswith("CPN60_PT")]
    pt = pt.set_index("proteotype_id")
    names = {"plastid_Cpn60beta": "plastid Cpn60β", "plastid_Cpn60alpha":
             "plastid Cpn60α", "mitochondrial_CPN60": "mitochondrial CPN60",
             "mitochondrial_groupI": "mitochondrial CPN60"}
    panels = []
    for cls in ("plastid_Cpn60beta", "plastid_Cpn60alpha", "mitochondrial_CPN60",
                "mitochondrial_groupI"):
        members = [p for p in ids if classes.get(p) == cls]
        if not members:
            continue

        def key(p):
            r = pt.loc[p]
            novel = str(r["novel_vs_reference"]).lower() in ("true", "1")
            return (-int(r["n_cultivars"]), novel, -int(r["length"]), p)
        members = sorted(members, key=key)
        ref = members[0]
        labels = {}
        for p in members:
            r = pt.loc[p]
            star = "*" if str(r["novel_vs_reference"]).lower() in ("true", "1") else " "
            labels[p] = f"{p.replace('CPN60_', '')} ({int(r['n_cultivars']):>2}){star}"
        P = build_panel(aln, members, ref, labels)
        P["name"] = (f"{names[cls]}: {len(members)} proteoforms against the most prevalent "
                     f"member of the class ({ref})")
        seq = P["seq"][ref].replace("-", "")
        feats = []
        m = re.search(r"GDGTT", seq)
        if m:
            a, b = col_of(P, ref, m.start() + 1), col_of(P, ref, m.end())
            feats.append({"label": "GDGTT nucleotide-binding loop", "c0": a, "c1": b,
                          "colour": "#2a78d6"})
        P["features_tracks"] = [("motifs", feats)]
        panels.append(P)
    return panels


def main(argv=None):
    global COLS, LINES, FONT, MONO
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", required=True, help="pan_proteome_build_* directory")
    ap.add_argument("--outdir", default=None,
                    help="default <inputs>/final/supplementary/figures")
    ap.add_argument("--cols", type=int, default=COLS)
    ap.add_argument("--lines", type=int, default=LINES)
    ap.add_argument("--font-size", type=float, default=FONT)
    ap.add_argument("--mono", default=MONO)
    a = ap.parse_args(argv)
    COLS, LINES, FONT, MONO = a.cols, a.lines, a.font_size, a.mono
    matplotlib.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
    outdir = a.outdir or os.path.join(a.inputs, "final", "supplementary", "figures")
    os.makedirs(outdir, exist_ok=True)

    rca_aln, cpn_aln = _find(a.inputs, "RCA.aln.fasta"), _find(a.inputs, "CPN60.aln.fasta")
    iso_p, ref_p = (_find(a.inputs, "rca_isoform_assignment.csv"),
                    _find(a.inputs, "rca_isoform_reference.json"))
    if rca_aln and iso_p and ref_p:
        P = rca_panel(read_aln(rca_aln), pd.read_csv(iso_p), json.load(open(ref_p)))
        n = render([P], "Residue-level alignment of the RCA proteoforms",
                   os.path.join(outdir, "RCA_alignment.pdf"),
                   os.path.join(outdir, "RCA_alignment_p1.png"),
                   "MAFFT alignment used for the RCA phylogeny, outgroups removed. Labels: "
                   "proteoform, established isoform (Supplementary Table of RCA isoform "
                   "assignments), v = variant gene model.")
        print(f"[alignment_views] RCA: {len(P['rows'])} rows, {len(P['num'])} columns, "
              f"{n} pages")
    else:
        print("[alignment_views] WARNING RCA view skipped: needs RCA.aln.fasta and the "
              "rca_isoforms.py outputs")
    if cpn_aln:
        cls = {}
        dc = _find(a.inputs, "domain_classification.csv")
        if dc:
            d = pd.read_csv(dc)
            cls = dict(zip(d["proteotype_id"], d["hmm_class"]))
        else:
            c2 = _find(a.inputs, "cpn60_classes.csv")
            if c2:
                d = pd.read_csv(c2)
                cls = dict(zip(d["proteotype_id"], d["cpn60_class"]))
        pt = pd.read_csv(_find(a.inputs, "proteotypes.csv"))
        panels = cpn60_panels(read_aln(cpn_aln), cls, pt)
        n = render(panels, "Residue-level alignment of the Cpn60 proteoforms by subunit class",
                   os.path.join(outdir, "CPN60_alignment.pdf"),
                   os.path.join(outdir, "CPN60_alignment_p1.png"),
                   "MAFFT alignment used for the Cpn60 phylogeny, outgroup removed, shown "
                   "per subunit class (HMM classification) with columns that are gaps in "
                   "every member of the class omitted.")
        print(f"[alignment_views] Cpn60: {sum(len(p['rows']) for p in panels)} rows in "
              f"{len(panels)} class panels, {n} pages")
    print(f"[alignment_views] outputs -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
