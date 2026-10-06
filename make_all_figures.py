#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_all_figures.py - regenerate every figure + data file of the wheat pan-proteome study."""

from __future__ import annotations
import math
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

plt.rcParams.update({
    "font.family": "DejaVu Sans", "axes.linewidth": 0.7,
    "axes.edgecolor": "#333333", "figure.dpi": 300, "savefig.dpi": 300,
})

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
except Exception:
    pass
BLUE, RED, GREY, ABSENT, AMBER = "#3B6FB6", "#C24B4B", "#B8BEC7", "#EFF1F3", "#E4A340"
FAMILIES = ("RCA", "CPN60")
FAM_DISP = {"RCA": "RCA", "CPN60": "Cpn60"}
CULT_LABELS = {
    "chinese_spring": "Chinese Spring", "arinalrfor": "ArinaLrFor", "jagger": "Jagger",
    "julius": "Julius", "lancer": "Lancer", "landmark": "Landmark", "mace": "Mace",
    "norin61": "Norin61", "stanley": "CDC Stanley", "mattis": "SY Mattis",
}
REFERENCE = "chinese_spring"
DEBUG_OVERLAP = False


def _audit(fig, ax, tag):
    """Count overlapping label pairs on `ax` (objective no-overlap QA)."""
    if not DEBUG_OVERLAP:
        return
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes = []
    for t in ax.texts:
        if not t.get_text().strip():
            continue
        try:
            boxes.append(t.get_window_extent(r))
        except Exception:
            pass
    bad = 0
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            dx = min(a.x1, b.x1) - max(a.x0, b.x0)
            dy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if dx > 1.5 and dy > 1.5:
                bad += 1
    print(f"[overlap-audit] {tag}: {len(boxes)} labels, {bad} overlapping pairs")


def clab(c):
    return CULT_LABELS.get(c, c.replace("_", " ").title())


def cultivar_columns(pav):
    return [c for c in pav.columns if c not in ("proteotype_id", "family")]


def _legend_below(fig, axes, handles, ncol=2, fontsize=10, pad=0.014):
    """Place a legend underneath everything that has already been drawn."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    axes = axes if isinstance(axes, (list, tuple)) else [axes]
    boxes = []
    for ax in axes:
        boxes.append(ax.get_window_extent(r))
        for t in list(ax.texts):
            if (t.get_text() or "").strip():
                boxes.append(t.get_window_extent(r))
        if getattr(ax, "axison", True):
            for t in ax.get_xticklabels():
                if (t.get_text() or "").strip():
                    boxes.append(t.get_window_extent(r))
            if (ax.get_xlabel() or "").strip():
                boxes.append(ax.xaxis.label.get_window_extent(r))
    y_px = min(b.y0 for b in boxes)
    y_fig = fig.transFigure.inverted().transform((0, y_px))[1]
    return fig.legend(handles=handles, loc="upper center",
                      bbox_to_anchor=(0.5, y_fig - pad), ncol=ncol,
                      frameon=False, fontsize=fontsize)


def _panel_letter(ax, letter, dx=-52, dy=16):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy),
                textcoords="offset points", fontsize=14, fontweight="bold",
                va="bottom", ha="left")


def plot_overview(pt, pav, cults, out_png):
    """Fig 1. (a) total vs novel proteotypes per family; (b) per-cultivar novel"""
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.2, 3.4))
    total = [int((pt.family == f).sum()) for f in FAMILIES]
    novel = [int(pt.loc[pt.family == f, "novel_vs_reference"].sum()) for f in FAMILIES]
    x = np.arange(len(FAMILIES)); w = 0.38
    axA.bar(x - w/2, total, w, label="Total proteotypes", color=BLUE, edgecolor="white")
    axA.bar(x + w/2, novel, w, label="Novel vs Chinese Spring", color=RED, edgecolor="white")
    for xi, (t, n) in enumerate(zip(total, novel)):
        axA.text(xi - w/2, t + 2, str(t), ha="center", va="bottom", fontsize=9)
        axA.text(xi + w/2, n + 2, f"{n} ({100*n/t:.0f}%)", ha="center", va="bottom", fontsize=8.5)
    axA.set_xticks(x); axA.set_xticklabels([FAM_DISP[f] for f in FAMILIES], fontsize=10)
    axA.set_ylabel("Non-redundant proteotypes", fontsize=10)
    axA.set_ylim(0, max(total) * 1.25); axA.tick_params(labelsize=9)
    axA.legend(frameon=False, fontsize=8.5, loc="upper left")
    axA.spines[["top", "right"]].set_visible(False)
    _panel_letter(axA, "a")

    def novel_counts(fam):
        ids = set(pt.loc[(pt.family == fam) & (pt.novel_vs_reference), "proteotype_id"])
        sub = pav[(pav.family == fam) & (pav.proteotype_id.isin(ids))]
        return {c: int(sub[c].sum()) for c in cults}
    cr, cc = novel_counts("RCA"), novel_counts("CPN60")
    order = sorted(cults, key=lambda c: -(cr[c] + cc[c]))
    xr = np.arange(len(order)); w = 0.4
    axB.bar(xr - w/2, [cr[c] for c in order], w, label="RCA", color=BLUE, edgecolor="white")
    axB.bar(xr + w/2, [cc[c] for c in order], w, label="Cpn60", color=AMBER, edgecolor="white")
    axB.set_xticks(xr)
    axB.set_xticklabels([clab(c) for c in order], rotation=45, ha="right", fontsize=8.5)
    axB.set_ylabel("Novel proteotypes carried", fontsize=10)
    axB.tick_params(labelsize=9)
    axB.legend(frameon=False, fontsize=9, loc="upper right")
    axB.spines[["top", "right"]].set_visible(False)
    _panel_letter(axB, "b")
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)


def _novel_matrix(pt, pav, family, cults, max_rows=None):
    novel = pt[(pt.family == family) & (pt.novel_vs_reference)]
    sub = pav[pav.proteotype_id.isin(set(novel.proteotype_id))].copy()
    sub["prev"] = sub[cults].sum(axis=1)
    sub = sub.sort_values(["prev", "proteotype_id"], ascending=[False, True])
    if max_rows and len(sub) > max_rows:
        sub = sub.head(max_rows)
    return sub["proteotype_id"].tolist(), sub[cults].to_numpy(int), int(len(novel))


def _draw_pav(ax, ids, M, cults, prefix, labels=None):
    ax.imshow(M, cmap=ListedColormap([ABSENT, BLUE]), aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(cults)))
    ax.set_xticklabels([clab(c) for c in cults], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(ids)))
    yl = [labels.get(i, i.replace(prefix, "")) for i in ids] if labels else [i.replace(prefix, "") for i in ids]
    ax.set_yticklabels(yl, fontsize=(5.6 if labels else 6.5))
    ax.tick_params(axis="y", pad=3)
    ax.set_xticks(np.arange(-.5, len(cults), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(ids), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="both", length=0)


def plot_novel_pav(pt, pav, cults, out_png, cpn60_rows=30, labels=None):
    """Fig 2. Presence/absence of NOVEL proteotypes: (a) RCA (all), (b) Cpn60"""
    rca_ids, rca_M, n_rca = _novel_matrix(pt, pav, "RCA", cults)
    cpn_ids, cpn_M, n_cpn = _novel_matrix(pt, pav, "CPN60", cults, max_rows=cpn60_rows)
    h = max(4.5, 0.17 * max(len(rca_ids), len(cpn_ids)) + 1.6)
    composite = bool(labels)
    fig = plt.figure(figsize=((10.4 if composite else 7.6), h))
    gs = fig.add_gridspec(1, 2, wspace=(0.85 if composite else 0.45))
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    _draw_pav(axA, rca_ids, rca_M, cults, "RCA_PT", labels=labels)
    axA.set_ylabel(f"{n_rca} novel RCA proteotypes", fontsize=9)
    _panel_letter(axA, "a")
    _draw_pav(axB, cpn_ids, cpn_M, cults, "CPN60_", labels=labels)
    axB.set_ylabel(f"{len(cpn_ids)} most prevalent of {n_cpn} novel Cpn60 proteotypes",
                   fontsize=9)
    _panel_letter(axB, "b")
    _legend_below(fig, [axA, axB],
                  [Patch(facecolor=BLUE, label="present"),
                   Patch(facecolor=ABSENT, edgecolor="#ccc", label="absent")],
                  ncol=2, fontsize=9)
    fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)


def plot_pav_heatmap(pt, pav, cults, family, out_png, labels=None):
    """Full presence/absence heatmap for one family. Input: pav_matrix.csv."""
    sub = pav[pav.family == family].reset_index(drop=True)
    ids = sub["proteotype_id"].tolist()
    novel = set(pt.loc[(pt.family == family) & (pt.novel_vs_reference), "proteotype_id"])
    M = sub[cults].to_numpy(int)
    fig, ax = plt.subplots(figsize=(6.4, max(4, 0.16 * len(ids) + 1.5)))
    ax.imshow(M, cmap=ListedColormap([ABSENT, BLUE]), aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(cults)))
    ax.set_xticklabels([clab(c) for c in cults], rotation=45, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(ids)))
    yl = [labels.get(i, i.replace(family + "_", "")) for i in ids] if labels else [i.replace(family + "_", "") for i in ids]
    ax.set_yticklabels(yl, fontsize=(4.8 if labels else 5.5))
    for tick, pid in zip(ax.get_yticklabels(), ids):
        if pid in novel:
            tick.set_color(RED)
    ax.set_xticks(np.arange(-.5, len(cults), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(ids), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="both", length=0)
    ax.legend(handles=[Patch(facecolor=BLUE, label="present"),
                       Patch(facecolor=ABSENT, edgecolor="#ccc", label="absent")],
              loc="lower center", bbox_to_anchor=(0.5, -0.12 - 2.0/len(ids)),
              ncol=2, frameon=False, fontsize=8)
    fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)


def plot_circos_pav(pt, pav, cults, family, out_png):
    """Higher-visibility radial (circos) presence/absence plot for one family."""
    sub = pav[pav.family == family].reset_index(drop=True)
    ids = sub["proteotype_id"].tolist()
    novel = set(pt.loc[(pt.family == family) & (pt.novel_vs_reference), "proteotype_id"])
    n = len(ids)
    if n == 0:
        return
    gap = math.radians(26)
    span = 2 * math.pi - gap
    start = math.pi / 2 + gap / 2
    theta = np.array([start + span * i / n for i in range(n)])
    width = span / n * 0.95
    r0 = 6.0
    ring_w = 0.92

    big = n > 60
    fig, ax = plt.subplots(figsize=(14, 14) if big else (10, 10),
                           subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("E"); ax.set_theta_direction(1)
    for i, cult in enumerate(cults):
        r = r0 + i * ring_w
        colors = []
        for pid in ids:
            present = int(sub.loc[sub.proteotype_id == pid, cult].iloc[0]) == 1
            colors.append((RED if pid in novel else BLUE) if present else "white")
        ax.bar(theta, [ring_w * 0.86] * n, width=width, bottom=r,
               color=colors, edgecolor="#d9d9d9", linewidth=0.3)
    r_out = r0 + len(cults) * ring_w
    lab_fs = 3.6 if big else 6.5
    for th, pid in zip(theta, ids):
        deg = math.degrees(th) % 360
        rot, ha = deg, "left"
        if 90 < deg < 270:
            rot, ha = deg - 180, "right"
        ax.text(th, r_out + 0.4, pid.replace(family + "_", ""),
                rotation=rot, rotation_mode="anchor", ha=ha, va="center", fontsize=lab_fs)
    for i, cult in enumerate(cults):
        ax.text(math.pi / 2, r0 + i * ring_w + ring_w * 0.42, clab(cult),
                ha="center", va="center", fontsize=8.5 if not big else 7.5,
                color="#222")
    ax.set_rmax(r_out + (3.2 if big else 2.4))
    ax.set_rticks([]); ax.set_xticks([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)
    ax.legend(handles=[Patch(facecolor=BLUE, label="present (in Chinese Spring)"),
                       Patch(facecolor=RED, label="present, novel vs Chinese Spring"),
                       Patch(facecolor="white", edgecolor="#bbb", label="absent")],
              loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=3,
              frameon=False, fontsize=10)
    _audit(fig, ax, f"circos_{family}")
    fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)


def plot_entropy(poly, family, out_png):
    """Readable positional-variability plot for one family (larger fonts, taller"""
    sub = poly[poly.family == family]
    if sub.empty:
        return
    mean_e = sub["entropy_bits"].mean()
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.bar(sub["aln_column"], sub["entropy_bits"], width=1.0, color=RED, alpha=0.9)
    ax.axhline(mean_e, color="#333", lw=1.0, ls="--")
    ax.text(sub["aln_column"].max(), mean_e, f"  mean = {mean_e:.2f} bits",
            va="center", ha="left", fontsize=10, color="#333")
    ax.set_xlabel("Alignment column (residue position)", fontsize=12)
    ax.set_ylabel("Shannon entropy (bits)", fontsize=12)
    ax.tick_params(labelsize=10)
    ax.set_xlim(0, sub["aln_column"].max() * 1.02)
    ax.margins(x=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)


def parse_newick(text):
    """Minimal Newick parser -> nested dict {name,length,support,children}."""
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]
    pos = 0

    def parse():
        nonlocal pos
        node = {"name": None, "length": 0.0, "support": None, "children": []}
        if pos < len(text) and text[pos] == "(":
            pos += 1
            while True:
                node["children"].append(parse())
                if text[pos] == ",":
                    pos += 1
                elif text[pos] == ")":
                    pos += 1
                    break
        start = pos
        while pos < len(text) and text[pos] not in ":,()":
            pos += 1
        label = text[start:pos].strip().strip("'\"")
        if node["children"]:
            node["support"] = label or None
        else:
            node["name"] = label
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",()":
                pos += 1
            try:
                node["length"] = float(text[start:pos])
            except ValueError:
                node["length"] = 0.0
        return node
    return parse()


def render_circular_tree(newick_text, out_png, title, novel_ids=None, prefix="",
                         label_fontsize=None, labels=None):
    """Circular (fan) phylogram rendered with matplotlib. Tips are coloured by"""
    novel_ids = novel_ids or set()
    root = parse_newick(newick_text)

    leaves = []
    def collect(n):
        if n["children"]:
            for c in n["children"]:
                collect(c)
        else:
            leaves.append(n)
    collect(root)
    n_leaves = len(leaves)
    if n_leaves < 3:
        return None

    gap = math.radians(18)
    span = 2 * math.pi - gap
    start = math.pi / 2 + gap / 2
    for i, lf in enumerate(leaves):
        lf["angle"] = start + span * i / (n_leaves - 1)

    maxr = [0.0]
    def assign(n, r0):
        r = r0 + max(0.0, n["length"])
        n["radius"] = r
        maxr[0] = max(maxr[0], r)
        if n["children"]:
            for c in n["children"]:
                assign(c, r)
            angs = [c["angle"] for c in n["children"]]
            n["angle"] = 0.5 * (min(angs) + max(angs))
    assign(root, 0.0)
    R = maxr[0] or 1.0

    big = n_leaves > 60
    fs = label_fontsize or (3.2 if big else 6.3)
    fig, ax = plt.subplots(figsize=(14, 14) if big else (9.5, 9.5))
    ax.set_aspect("equal"); ax.axis("off")
    col = "#444"

    def draw(n, r_par):
        a = n["angle"]
        ax.plot([r_par * math.cos(a), n["radius"] * math.cos(a)],
                [r_par * math.sin(a), n["radius"] * math.sin(a)],
                color=col, lw=0.6)
        if n["children"]:
            angs = [c["angle"] for c in n["children"]]
            arc = np.linspace(min(angs), max(angs), 60)
            ax.plot(n["radius"] * np.cos(arc), n["radius"] * np.sin(arc), color=col, lw=0.6)
            for c in n["children"]:
                draw(c, n["radius"])
        else:
            deg = math.degrees(a) % 360
            rot, ha = deg, "left"
            if 90 < deg < 270:
                rot, ha = deg - 180, "right"
            name = labels.get(n["name"], (n["name"] or "").replace(prefix, "")) if labels else (n["name"] or "").replace(prefix, "")
            is_novel = (n["name"] in novel_ids)
            ax.plot([n["radius"] * math.cos(a)], [n["radius"] * math.sin(a)],
                    marker="o", ms=2.2, color=(RED if is_novel else BLUE))
            r_lab = R * 1.06
            ax.plot([n["radius"] * math.cos(a), r_lab * math.cos(a)],
                    [n["radius"] * math.sin(a), r_lab * math.sin(a)],
                    color="#cccccc", lw=0.3, zorder=0)
            rx, ry = (r_lab + 0.015 * R) * math.cos(a), (r_lab + 0.015 * R) * math.sin(a)
            ax.text(rx, ry, name, rotation=rot, rotation_mode="anchor",
                    ha=ha, va="center", fontsize=fs,
                    color=(RED if is_novel else "#222"))
    draw(root, 0.0)
    lim = R * 1.55
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    if novel_ids:
        _legend_below(fig, ax,
                      [Patch(facecolor=BLUE, label="in Chinese Spring"),
                       Patch(facecolor=RED, label="novel vs Chinese Spring")],
                      ncol=2, fontsize=10)
    _audit(fig, ax, title)
    fig.savefig(out_png, bbox_inches="tight", dpi=300); plt.close(fig)
    return out_png


def pangenome_partition(pav, cults):
    """Classify each cluster as core / soft-core / shell / cloud. Returns"""
    M = pav[cults].to_numpy(bool)
    prev = M.sum(1)
    N = len(cults)
    soft = math.ceil(0.9 * N)
    part = {
        "core": int((prev == N).sum()),
        "soft-core": int(((prev >= soft) & (prev < N)).sum()),
        "shell": int(((prev >= 2) & (prev < soft)).sum()),
        "cloud": int((prev == 1).sum()),
    }
    return part, prev


def rarefaction(pav, cults, nperm=50, seed=0):
    """Pan- and core-proteome accumulation over random cultivar orderings."""
    M = pav[cults].to_numpy(bool)
    N = len(cults)
    rng = np.random.default_rng(seed)
    pan = np.zeros((nperm, N)); core = np.zeros((nperm, N))
    idx = np.arange(N)
    for p in range(nperm):
        rng.shuffle(idx)
        cum_any = np.zeros(M.shape[0], bool)
        cum_all = np.ones(M.shape[0], bool)
        for k in range(N):
            col = M[:, idx[k]]
            cum_any |= col; cum_all &= col
            pan[p, k] = cum_any.sum(); core[p, k] = cum_all.sum()
    return pan, core


def plot_pangenome_overview(pav, cults, out_png,
                            title="Whole-proteome pan-genome (Chinese Spring reference)"):
    """MAIN genome-wide figure: (a) core/soft-core/shell/cloud partition;"""
    part, prev = pangenome_partition(pav, cults)
    pan, core = rarefaction(pav, cults)
    N = len(cults)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    cats = ["core", "soft-core", "shell", "cloud"]
    vals = [part[c] for c in cats]; tot = max(1, sum(vals))
    colors = [BLUE, "#6C9BD1", AMBER, RED]
    axA.bar(range(len(cats)), vals, color=colors, edgecolor="white")
    for i, v in enumerate(vals):
        axA.text(i, v + tot * 0.012, f"{v}\n({100*v/tot:.0f}%)", ha="center", va="bottom", fontsize=8)
    axA.set_xticks(range(len(cats)))
    axA.set_xticklabels([f"Core\n(all {N})", "Soft-core\n(\u226590%)", "Shell", "Cloud\n(private)"], fontsize=8.5)
    axA.set_ylabel("Proteotype clusters", fontsize=10)
    axA.set_ylim(0, max(vals) * 1.22); axA.tick_params(labelsize=9)
    axA.spines[["top", "right"]].set_visible(False); _panel_letter(axA, "a")
    ks = np.arange(1, N + 1)
    for arr, lab, c in [(pan, "Pan-proteome", BLUE), (core, "Core proteome", RED)]:
        axB.plot(ks, arr.mean(0), "-o", color=c, ms=3.5, label=lab)
        axB.fill_between(ks, arr.min(0), arr.max(0), color=c, alpha=0.15)
    axB.set_xlabel("Cultivars sampled", fontsize=10)
    axB.set_ylabel("Proteotype clusters", fontsize=10)
    axB.set_xticks(ks); axB.tick_params(labelsize=9)
    axB.legend(frameon=False, fontsize=9); axB.spines[["top", "right"]].set_visible(False)
    _panel_letter(axB, "b")
    fig.tight_layout(); fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)
    return part


def plot_prevalence_hist(pav, cults, out_png):
    """SUPP genome-wide figure: prevalence spectrum (U-shape). Input: genome_wide_pav.csv."""
    _, prev = pangenome_partition(pav, cults)
    N = len(cults)
    fig, ax = plt.subplots(figsize=(7, 3.4))
    ax.hist(prev, bins=np.arange(0.5, N + 1.5, 1), color=BLUE, edgecolor="white")
    ax.set_xlabel("Number of cultivars in which a proteotype is present", fontsize=11)
    ax.set_ylabel("Number of proteotype clusters", fontsize=11)
    ax.set_xticks(range(1, N + 1)); ax.tick_params(labelsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)


def plot_genomewide_per_cultivar(pav, cults, out_png):
    """SUPP genome-wide figure: per-cultivar core vs accessory composition."""
    M = pav[cults].to_numpy(bool); prev = M.sum(1); N = len(cults)
    core_mask = prev == N
    core_counts = [int((M[:, j] & core_mask).sum()) for j in range(N)]
    acc_counts = [int((M[:, j] & ~core_mask).sum()) for j in range(N)]
    order = sorted(range(N), key=lambda j: -(core_counts[j] + acc_counts[j]))
    x = np.arange(N)
    cc = [core_counts[j] for j in order]; ac = [acc_counts[j] for j in order]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.bar(x, cc, color=BLUE, label="core", edgecolor="white")
    ax.bar(x, ac, bottom=cc, color=AMBER, label="accessory", edgecolor="white")
    ax.set_xticks(x); ax.set_xticklabels([clab(cults[j]) for j in order], rotation=45, ha="right", fontsize=8.5)
    ax.set_ylabel("Proteotype clusters", fontsize=10); ax.tick_params(labelsize=9)
    ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)


def parse_uniprot_id(qseqid):
    """UniProt accession + entry name from a DIAMOND qseqid such as"""
    parts = str(qseqid).split("|")
    if len(parts) >= 3:
        return parts[1], parts[2]
    return str(qseqid), ""


def build_proteotype_naming(inputs_dir, out_csv=None):
    """Nomenclature/annotation table linking every proteotype to its closest"""
    import os
    pt = pd.read_csv(os.path.join(inputs_dir, "proteotypes.csv"))
    mem = pd.read_csv(os.path.join(inputs_dir, "members.csv"))
    hits = pd.read_csv(os.path.join(inputs_dir, "diamond_hits_all.csv"))
    hits = hits.sort_values(["evalue", "bitscore"], ascending=[True, False])
    best = hits.groupby(["cultivar", "sseqid"]).head(1).copy()
    parsed = [parse_uniprot_id(q) for q in best["qseqid"]]
    best["uacc"] = [a for a, _ in parsed]
    best["uentry"] = [e for _, e in parsed]
    key2u = {(r.cultivar, r.sseqid): (r.uacc, r.uentry, r.evalue) for r in best.itertuples()}
    rows = []
    for m in mem.itertuples():
        u = key2u.get((m.cultivar, m.source_id), ("", "", 1.0))
        rows.append((m.proteotype_id, m.cultivar, m.source_id, u[0], u[1], u[2],
                     getattr(m, "evalue_to_ref", 1.0)))
    md = pd.DataFrame(rows, columns=["proteotype_id", "cultivar", "source_id",
                                     "uacc", "uentry", "hit_evalue", "evalue_to_ref"])

    def agg(g):
        if (g.cultivar == "chinese_spring").any():
            rep = g[g.cultivar == "chinese_spring"].sort_values(
                "source_id", kind="mergesort").iloc[0]
        else:
            rep = g.sort_values(
                ["evalue_to_ref", "cultivar", "source_id"], kind="mergesort").iloc[0]
        ua = g[g.uacc != ""].sort_values(
            ["hit_evalue", "uacc"], kind="mergesort")
        uacc = ua["uacc"].iloc[0] if len(ua) else ""
        uentry = ua["uentry"].iloc[0] if len(ua) else ""
        return pd.Series({"representative_cultivar": rep.cultivar,
                          "representative_source_id": rep.source_id,
                          "uniprot_id": uacc, "uniprot_entry": uentry})
    _agg_cols = ["cultivar", "source_id", "uacc", "uentry", "hit_evalue", "evalue_to_ref"]
    rep = md.groupby("proteotype_id")[_agg_cols].apply(agg).reset_index()
    nm = pt.merge(rep, on="proteotype_id", how="left")

    def comp(r):
        acc = r.uniprot_id if isinstance(r.uniprot_id, str) and r.uniprot_id else "NA"
        return f"{r.proteotype_id}__{r.representative_cultivar}__{acc}"
    nm["composite_id"] = nm.apply(comp, axis=1)
    keep = ["proteotype_id", "family", "uniprot_id", "uniprot_entry",
            "representative_cultivar", "representative_source_id", "n_cultivars",
            "novel_vs_reference", "composite_id", "cultivars", "source_ids"]
    nm = nm[[c for c in keep if c in nm.columns]]
    if out_csv:
        nm.to_csv(out_csv, index=False)
    return nm


def composite_label_map(naming_df, style="compact"):
    """proteotype_id -> display label. compact -> 'RCA_PT001|A0A..|Chinese Spring'."""
    out = {}
    for r in naming_df.itertuples():
        cult = CULT_LABELS.get(r.representative_cultivar, r.representative_cultivar)
        acc = r.uniprot_id if isinstance(r.uniprot_id, str) and r.uniprot_id else "NA"
        out[r.proteotype_id] = (r.composite_id if style == "full"
                                else f"{r.proteotype_id}|{acc}|{cult}")
    return out


def write_annotated_fasta(sequences_fa, naming_df, out_fa):
    """Re-header a proteotype FASTA (record ids = proteotype_id) with composite"""
    info = {r.proteotype_id: r for r in naming_df.itertuples()}

    def header(pid):
        r = info.get(pid)
        if r is None:
            return pid, ""
        desc = (f"family={r.family} uniprot={r.uniprot_id} uniprot_entry={r.uniprot_entry} "
                f"rep_cultivar={r.representative_cultivar} rep_source_id={r.representative_source_id} "
                f"n_cultivars={r.n_cultivars} novel={r.novel_vs_reference}")
        return r.composite_id, desc

    n = 0
    with open(sequences_fa) as fi, open(out_fa, "w") as fo:
        cur, seq = None, []

        def flush():
            nonlocal n
            if cur is None:
                return
            pid = cur.split()[0]
            cid, desc = header(pid)
            fo.write(f">{cid} {desc}\n")
            s = "".join(seq)
            for i in range(0, len(s), 60):
                fo.write(s[i:i + 60] + "\n")
            n += 1
        for line in fi:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush(); cur = line[1:]; seq = []
            else:
                seq.append(line)
        flush()
    return out_fa, n


import argparse
import os
import shutil

FAMILY_PREFIX = {"RCA": "RCA_", "CPN60": "CPN60_"}


def _find_input(inputs_dir, filename):
    """Locate a pipeline output. Bulky tables live in <result>/supplementary/, the rest"""
    for cand in (os.path.join(inputs_dir, filename),
                 os.path.join(inputs_dir, "supplementary", filename)):
        if os.path.exists(cand):
            return cand
    return None


def _read_nwk(path):
    return open(path).read() if os.path.exists(path) else None


def _write_supp_index(suppdir, made):
    figs = sorted(m for m in made if m.endswith(".png"))
    tabs = sorted(m for m in made if m.endswith(".csv"))
    lines = ["# Supplementary figures and data", "",
             "Generated by make_all_figures.py.", ""]
    if figs:
        lines.append("## Supplementary figures")
        lines += [f"- **{m}**" for m in figs] + [""]
    if tabs:
        lines.append("## Supplementary data tables")
        lines += [f"- **{m}**" for m in tabs] + [""]
    open(os.path.join(suppdir, "SUPPLEMENTARY_INDEX.md"), "w").write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", required=True,
                    help="pan_proteome_build_* dir with CSVs (+ optional *.nwk, genome_wide_pav.csv).")
    ap.add_argument("--outdir", default="./figures", help="Main figures (default ./figures).")
    ap.add_argument("--suppdir", default=None,
                    help="Supplementary figures + data (default <outdir>/supplementary).")
    ap.add_argument("--cpn60-rows", type=int, default=30,
                    help="Cpn60 novel proteotypes shown in Figure 2b (default 30).")
    ap.add_argument("--sequences", default=None,
                    help="Proteotype FASTA (record ids = proteotype ids) -> annotated FASTA.")
    ap.add_argument("--label", choices=["id", "composite"], default="id",
                    help="Labels on PAV heatmaps + tree tips: short id (default) or composite.")
    a = ap.parse_args()
    outdir = a.outdir
    suppdir = a.suppdir or os.path.join(outdir, "supplementary")
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(suppdir, exist_ok=True)

    pt = pd.read_csv(os.path.join(a.inputs, "proteotypes.csv"))
    pt["novel_vs_reference"] = pt["novel_vs_reference"].astype(str).str.lower().isin(["true", "1"])
    pav = pd.read_csv(os.path.join(a.inputs, "pav_matrix.csv"))
    poly = pd.read_csv(os.path.join(a.inputs, "polymorphic_sites.csv"))
    cults = cultivar_columns(pav)
    made_supp = []

    labels = None
    naming_src = a.inputs
    if not all(os.path.exists(os.path.join(naming_src, f))
               for f in ("proteotypes.csv", "members.csv", "diamond_hits_all.csv")):
        alt = os.path.join(a.inputs, "supplementary")
        if all(os.path.exists(os.path.join(alt, f))
               for f in ("proteotypes.csv", "members.csv", "diamond_hits_all.csv")):
            naming_src = alt
    have_naming = all(os.path.exists(os.path.join(naming_src, f))
                      for f in ("proteotypes.csv", "members.csv", "diamond_hits_all.csv"))
    if have_naming:
        nm = build_proteotype_naming(naming_src, os.path.join(outdir, "proteotype_naming.csv"))
        nm.to_csv(os.path.join(suppdir, "supp_proteotype_naming.csv"), index=False)
        made_supp.append("supp_proteotype_naming.csv")
        if a.label == "composite":
            labels = composite_label_map(nm)
        if a.sequences and os.path.exists(a.sequences):
            out_fa = os.path.join(outdir, "PanProteome_catalog.annotated.faa")
            _, ncount = write_annotated_fasta(a.sequences, nm, out_fa)
            print(f"[naming] annotated FASTA: {ncount} records -> {out_fa}")
        elif a.sequences:
            print(f"[naming] --sequences {a.sequences} not found; annotated FASTA skipped.")
        print(f"[naming] proteotype_naming.csv written ({len(nm)} proteotypes).")
    else:
        print("[naming] members.csv / diamond_hits_all.csv not found -> naming table + "
              "annotated FASTA skipped; figures use short ids.")

    plot_overview(pt, pav, cults, os.path.join(outdir, "Figure1_overview.png"))
    plot_novel_pav(pt, pav, cults, os.path.join(outdir, "Figure2_novel_pav.png"),
                   cpn60_rows=a.cpn60_rows, labels=labels)
    for fam, fname in [("RCA", "Figure3_RCA_tree_circular.png"),
                       ("CPN60", "Figure4_CPN60_tree_circular.png")]:
        nwk = _read_nwk(os.path.join(a.inputs, f"{fam}.nwk"))
        if nwk:
            novel = set(pt.loc[(pt.family == fam) & pt.novel_vs_reference, "proteotype_id"])
            done = False
            try:
                import tree_tools
                troot = tree_tools.parse_newick(nwk)
                summ = tree_tools.support_summary(troot)
                if summ.get("n"):
                    print(f"[tree] {fam}: {summ['n']} internal nodes, median support "
                          f"{summ['median']:.3f}, {summ['ge_0.80']} at >= 0.80")
                thr = float(os.environ.get("PANPROTEOME_COLLAPSE", "0.8"))
                n_col = tree_tools.collapse_unsupported(troot, thr)
                troot = tree_tools.midpoint_root(troot)
                with open(os.path.join(outdir, f"{fam}.supported.nwk"), "w") as fh:
                    fh.write(tree_tools.to_newick(troot) + ";\n")
                tree_tools.render_circular(
                    troot, os.path.join(outdir, fname), novel_ids=novel,
                    prefix=FAMILY_PREFIX[fam],
                    title_note=(f"Node labels: SH-like local support; nodes < {thr:g} "
                                f"collapsed ({n_col}); midpoint-rooted"))
                done = True
            except Exception as e:
                print(f"[tree] tree_tools unavailable/failed for {fam} ({e}); "
                      f"using the plain renderer.")
            if not done:
                render_circular_tree(nwk, os.path.join(outdir, fname),
                                     f"{FAM_DISP[fam]} proteotype tree (circular)",
                                     novel_ids=novel, prefix=FAMILY_PREFIX[fam], labels=labels)
        else:
            print(f"[main] {fam}.nwk not found -> circular tree skipped (the pipeline writes it).")
    plot_entropy(poly, "RCA", os.path.join(outdir, "Figure6_RCA_variability.png"))
    plot_entropy(poly, "CPN60", os.path.join(outdir, "Figure7_CPN60_variability.png"))

    gw_path = _find_input(a.inputs, "genome_wide_pav.csv")
    if gw_path:
        gw = pd.read_csv(gw_path)
        gwc = [c for c in gw.columns if c not in ("cluster_id", "family", "novel_vs_reference")]
        part = plot_pangenome_overview(gw, gwc, os.path.join(outdir, "Figure5_pangenome_overview.png"))
        plot_prevalence_hist(gw, gwc, os.path.join(suppdir, "FigureS5_prevalence_spectrum.png"))
        plot_genomewide_per_cultivar(gw, gwc, os.path.join(suppdir, "FigureS6_per_cultivar_composition.png"))
        made_supp += ["FigureS5_prevalence_spectrum.png", "FigureS6_per_cultivar_composition.png"]
        _, prev = pangenome_partition(gw, gwc)
        pd.DataFrame({"category": list(part), "n_clusters": list(part.values())}).to_csv(
            os.path.join(suppdir, "supp_pangenome_partition.csv"), index=False)
        gw.assign(n_cultivars=prev).to_csv(
            os.path.join(suppdir, "supp_genome_wide_pav.csv"), index=False)
        made_supp += ["supp_genome_wide_pav.csv", "supp_pangenome_partition.csv"]
    else:
        print("[genome-wide] genome_wide_pav.csv not found -> Figure 5 and S5/S6 skipped.")
        print("              (the pipeline writes it to <result>/supplementary/ when")
        print("               whole_proteome.enabled: true in config.yaml)")

    hm = {"RCA": "FigureS1_RCA_PAV_heatmap.png", "CPN60": "FigureS2_CPN60_PAV_heatmap.png"}
    ci = {"RCA": "circos_RCA_PAV.png", "CPN60": "circos_CPN60_PAV.png"}
    for fam in FAMILIES:
        plot_pav_heatmap(pt, pav, cults, fam, os.path.join(suppdir, hm[fam]), labels=labels); made_supp.append(hm[fam])
        if os.environ.get("PANPROTEOME_DRAW_CIRCOS"):
            plot_circos_pav(pt, pav, cults, fam, os.path.join(suppdir, ci[fam]))
            made_supp.append(ci[fam])

    for src_name, dst in [("proteotypes.csv", "supp_proteotypes.csv"),
                          ("pav_matrix.csv", "supp_pav_matrix.csv"),
                          ("physicochemical_features.csv", "supp_physicochemical_features.csv"),
                          ("polymorphic_sites.csv", "supp_polymorphic_sites.csv")]:
        src = _find_input(a.inputs, src_name)
        dst_path = os.path.join(suppdir, dst)
        if src and os.path.abspath(src) != os.path.abspath(dst_path):
            shutil.copy(src, dst_path); made_supp.append(dst)

    _write_supp_index(suppdir, made_supp)
    print("MAIN outputs      ->", os.path.abspath(outdir))
    print("SUPPLEMENTARY dir ->", os.path.abspath(suppdir))


if __name__ == "__main__":
    main()
