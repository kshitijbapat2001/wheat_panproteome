#!/usr/bin/env python3
import argparse
import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Wedge

MM = 1 / 25.4
OI = {"blue": "#0072B2", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
      "yellow": "#F0E442", "vermilion": "#D55E00", "purple": "#CC79A7",
      "grey": "#999999", "light": "#C8C8C8", "absent": "#EFEFEF"}
CULT = {"chinese_spring": "Chinese Spring", "arinalrfor": "ArinaLrFor", "jagger": "Jagger",
        "julius": "Julius", "lancer": "Lancer", "landmark": "Landmark", "mace": "Mace",
        "norin61": "Norin61", "stanley": "CDC Stanley", "mattis": "SY Mattis"}
CULT_ORDER = ["chinese_spring", "arinalrfor", "jagger", "julius", "lancer", "landmark",
              "mace", "mattis", "norin61", "stanley"]

matplotlib.rcParams.update({
    "font.family": "Arial", "font.size": 6, "axes.labelsize": 6.5, "axes.titlesize": 6.5,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})


def find(root, name):
    hits = [os.path.join(d, name) for d, _, files in os.walk(root) if name in files
            and "tables_formatted" not in d]
    if not hits:
        raise FileNotFoundError(name)
    return sorted(hits, key=lambda h: (h.count(os.sep), "final" in h, h))[0]


def truthy(s):
    return s.astype(str).str.lower().isin(["true", "1"])


def label(fig, x_mm, y_mm, letter, H):
    fig.text(x_mm * MM / fig.get_figwidth(), 1 - y_mm * MM / fig.get_figheight(), letter,
             fontsize=8, fontweight="bold", va="top", ha="left")


def axes_mm(fig, x, y, w, h, H, W=180):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def save(fig, out, name):
    fig.savefig(os.path.join(out, name + ".pdf"))
    fig.savefig(os.path.join(out, name + ".png"), dpi=600)
    plt.close(fig)


def thousands(ax, axis="y"):
    fmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}")
    (ax.yaxis if axis == "y" else ax.xaxis).set_major_formatter(fmt)


def fig1(pt, pav, out):
    W, H = 180, 62
    fig = plt.figure(figsize=(W * MM, H * MM))
    fams = [("RCA", "RCA"), ("CPN60", "Cpn60")]
    tot = [int((pt.family == f).sum()) for f, _ in fams]
    nov = [int((pt.family == f).loc[pt.novel].sum()) for f, _ in fams]
    ax = axes_mm(fig, 14, 6, 52, 44, H)
    x = np.arange(2)
    w = 0.36
    ax.bar(x - w / 2, tot, w, color=OI["light"], label="All proteoforms", linewidth=0)
    ax.bar(x + w / 2, nov, w, color=OI["vermilion"], label="Reference-absent", linewidth=0)
    for i in range(2):
        ax.text(x[i] - w / 2, tot[i] + 2, f"{tot[i]}", ha="center", va="bottom", fontsize=5.5)
        ax.text(x[i] + w / 2, nov[i] + 2, f"{nov[i]} ({100 * nov[i] / tot[i]:.0f}%)",
                ha="center", va="bottom", fontsize=5.5)
    ax.set_xticks(x)
    ax.set_xticklabels([d for _, d in fams])
    ax.set_ylabel("Proteoforms (n)")
    ax.set_ylim(0, 165)
    ax.legend(loc="upper left", handlelength=1, handleheight=0.8)
    label(fig, 1, 2, "a", H)

    cults = CULT_ORDER
    counts = {}
    for f, _ in fams:
        ids = set(pt.loc[(pt.family == f) & pt.novel, "proteotype_id"])
        sub = pav[(pav.family == f) & pav.proteotype_id.isin(ids)]
        counts[f] = {c: int(sub[c].astype(int).sum()) for c in cults}
    order = sorted(cults, key=lambda c: (-(counts["RCA"][c] + counts["CPN60"][c]),
                                         -counts["RCA"][c], CULT[c]))
    ax = axes_mm(fig, 88, 6, 88, 36, H)
    xr = np.arange(len(order))
    ax.bar(xr - 0.2, [counts["RCA"][c] for c in order], 0.4, color=OI["blue"], label="RCA",
           linewidth=0)
    ax.bar(xr + 0.2, [counts["CPN60"][c] for c in order], 0.4, color=OI["orange"],
           label="Cpn60", linewidth=0)
    ax.set_xticks(xr)
    ax.set_xticklabels([CULT[c] for c in order], rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel("Reference-absent\nproteoforms (n)")
    ax.legend(loc="upper right", handlelength=1, handleheight=0.8)
    label(fig, 76, 2, "b", H)
    save(fig, out, "Fig1")
    rows_a = pd.DataFrame({"family": [d for _, d in fams], "proteoforms": tot,
                           "reference_absent": nov})
    rows_b = pd.DataFrame({"cultivar": [CULT[c] for c in order],
                           "RCA_reference_absent": [counts["RCA"][c] for c in order],
                           "Cpn60_reference_absent": [counts["CPN60"][c] for c in order]})
    return {"Fig1a": rows_a, "Fig1b": rows_b}


class Node:
    __slots__ = ("name", "length", "support", "children", "parent", "r", "theta")

    def __init__(self):
        self.name, self.length, self.support = "", 0.0, None
        self.children, self.parent, self.r, self.theta = [], None, 0.0, 0.0


def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                c = node()
                c.parent = n
                n.children.append(c)
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        m = re.match(r"[^:,()]*", s[pos:])
        tok = m.group(0)
        pos += len(tok)
        if n.children:
            try:
                n.support = float(tok) if tok else None
            except ValueError:
                n.name = tok
        else:
            n.name = tok
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m = re.match(r"[-0-9.eE]+", s[pos:])
            n.length = float(m.group(0))
            pos += len(m.group(0))
        return n
    return node()


def leaves(n):
    return [n] if not n.children else [x for c in n.children for x in leaves(c)]


def nodes(n):
    return [n] + [x for c in n.children for x in nodes(c)]


def prune(root, names):
    for leaf in [l for l in leaves(root) if l.name in names]:
        par = leaf.parent
        par.children.remove(leaf)
        while par is not None and len(par.children) == 1:
            child = par.children[0]
            if par.parent is None:
                child.parent, child.length, root = None, 0.0, child
            else:
                child.length += par.length
                child.parent = par.parent
                par.parent.children[par.parent.children.index(par)] = child
            par = child.parent
    return root


def layout(root, gap_deg=14.0):
    lv = leaves(root)
    step = (360.0 - gap_deg) / len(lv)
    for i, leaf in enumerate(lv):
        leaf.theta = 90 - gap_deg / 2 - i * step

    def depth(n, d):
        n.r = d
        for c in n.children:
            depth(c, d + max(c.length, 0.0))
    depth(root, 0.0)

    def ang(n):
        if n.children:
            for c in n.children:
                ang(c)
            n.theta = (n.children[0].theta + n.children[-1].theta) / 2
    ang(root)
    return lv, step


def draw_tree(ax, root, tip_status, tip_class, class_colors, label_fn, R, label_fs=5,
              scale=0.1, breaks=None):
    broken = []
    for n, show in (breaks or {}).items():
        broken.append((n, n.length))
        n.length = show
    lv, step = layout(root)
    maxr = max(l.r for l in lv)
    k = R / maxr

    def xy(r, th):
        return r * math.cos(math.radians(th)), r * math.sin(math.radians(th))

    def walk(n):
        for c in n.children:
            x0, y0 = xy(n.r * k, c.theta)
            x1, y1 = xy(c.r * k, c.theta)
            ax.plot([x0, x1], [y0, y1], color="black", lw=0.45, solid_capstyle="butt")
            if any(c is b for b, _ in broken):
                xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
                t = math.radians(c.theta)
                for d in (-0.012, 0.012):
                    cx, cy = xm + d * R * math.cos(t), ym + d * R * math.sin(t)
                    u = t + math.radians(60)
                    ax.plot([cx - 0.03 * R * math.cos(u), cx + 0.03 * R * math.cos(u)],
                            [cy - 0.03 * R * math.sin(u), cy + 0.03 * R * math.sin(u)],
                            color="black", lw=0.5, zorder=6)
                ax.plot([xm - 0.012 * R * math.cos(t), xm + 0.012 * R * math.cos(t)],
                        [ym - 0.012 * R * math.sin(t), ym + 0.012 * R * math.sin(t)],
                        color="white", lw=1.2, zorder=5)
            walk(c)
        if n.children:
            t0, t1 = n.children[-1].theta, n.children[0].theta
            ts = np.linspace(t0, t1, max(2, int(abs(t1 - t0) * 2)))
            ax.plot(n.r * k * np.cos(np.radians(ts)), n.r * k * np.sin(np.radians(ts)),
                    color="black", lw=0.45, solid_capstyle="butt")
            if n.support is not None and n.parent is not None:
                x, y = xy(n.r * k, n.theta)
                ax.plot(x, y, "o", ms=1.9, mew=0,
                        color="black" if n.support >= 0.95 else OI["grey"], zorder=4)
    walk(root)
    r_lab = R * 1.07
    r_ring = R * 1.035
    for l in lv:
        x, y = xy(l.r * k, l.theta)
        st = tip_status.get(l.name)
        if st is None:
            ax.plot(x, y, "o", ms=2.2, mfc="white", mec="black", mew=0.4, zorder=5)
        else:
            ax.plot(x, y, "o", ms=2.2, mew=0, zorder=5,
                    color=OI["vermilion"] if st else OI["blue"])
        xl, yl = xy(r_ring - R * 0.012, l.theta)
        ax.plot([x, xl], [y, yl], color="#BDBDBD", lw=0.25, zorder=1)
        cls = tip_class.get(l.name)
        if cls:
            ax.add_patch(Wedge((0, 0), r_ring + R * 0.022, l.theta - step / 2,
                               l.theta + step / 2, width=R * 0.03, color=class_colors[cls],
                               lw=0))
        xt, yt = xy(r_lab, l.theta)
        th = l.theta % 360
        flip = 90 < th < 270
        ax.text(xt, yt, label_fn(l.name), fontsize=label_fs, ha="right" if flip else "left",
                va="center", rotation=th + 180 if flip else th, rotation_mode="anchor")
    x0 = -R * 1.25
    y0 = -R * 1.28
    ax.plot([x0, x0 + scale * k], [y0, y0], color="black", lw=0.6)
    ax.text(x0 + scale * k / 2, y0 - R * 0.035, f"{scale}", ha="center", va="top",
            fontsize=5.5)
    lim = R * 1.42
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    for n, true_len in broken:
        n.length = true_len


def fig2(build, pt, iso, dom, out):
    rca = prune(parse_newick(open(find(build, "RCA.supported.nwk")).read()),
                {"OUTGROUP_Q9SB58"})
    cpn = parse_newick(open(find(build, "CPN60.supported.nwk")).read())
    status = dict(zip(pt.proteotype_id, pt.novel))
    rca_cls = {r.proteotype_id: r.isoform for r in iso.itertuples()}
    rca_colors = {"RCA1β": OI["green"], "RCA2α": OI["orange"],
                  "RCA2β": OI["sky"], "RCA-like": OI["purple"]}
    cname = {"plastid_Cpn60alpha": "Plastid Cpn60α", "plastid_Cpn60beta":
             "Plastid Cpn60β", "mitochondrial_CPN60": "Mitochondrial CPN60"}
    cpn_cls = {r.proteotype_id: cname.get(r.hmm_class) for r in dom.itertuples()
               if r.proteotype_id.startswith("CPN60_")}
    cpn_colors = {"Plastid Cpn60α": OI["orange"], "Plastid Cpn60β": OI["green"],
                  "Mitochondrial CPN60": OI["purple"]}

    def lab(n):
        return n.replace("OUTGROUP_", "").replace("RCA_", "").replace("CPN60_", "")
    W, H = 180, 112
    fig = plt.figure(figsize=(W * MM, H * MM))
    rca_like = set(iso.loc[iso.isoform == "RCA-like", "proteotype_id"])
    stem = [n for n in nodes(rca) if n.children and {l.name for l in leaves(n)} == rca_like]
    axa = axes_mm(fig, 0, 4, 72, 72, H)
    draw_tree(axa, rca, status, rca_cls, rca_colors, lab, R=1.0, label_fs=5.2,
              breaks={stem[0]: 0.2} if stem else None)
    axb = axes_mm(fig, 70, 0, 110, 110, H)
    draw_tree(axb, cpn, status, cpn_cls, cpn_colors, lab, R=1.0, label_fs=5.0)
    label(fig, 1, 2, "a", H)
    label(fig, 72, 2, "b", H)
    common = [Line2D([], [], marker="o", ls="", ms=3, mew=0, color=OI["blue"],
                     label="Present in Chinese Spring"),
              Line2D([], [], marker="o", ls="", ms=3, mew=0, color=OI["vermilion"],
                     label="Reference-absent"),
              Line2D([], [], marker="o", ls="", ms=3, mfc="white", mec="black", mew=0.5,
                     label="Outgroup"),
              Line2D([], [], marker="o", ls="", ms=2.5, mew=0, color="black",
                     label="Local support ≥0.95"),
              Line2D([], [], marker="o", ls="", ms=2.5, mew=0, color=OI["grey"],
                     label="Local support 0.80–0.95")]
    fig.legend(handles=common, loc="upper left", bbox_to_anchor=(2 / W, 1 - 80 / H),
               ncol=1, handletextpad=0.3, borderaxespad=0)
    fig.legend(handles=[Patch(color=c, label=k) for k, c in rca_colors.items()],
               loc="upper left", bbox_to_anchor=(38 / W, 1 - 80 / H), ncol=1,
               handlelength=1, handleheight=0.7, handletextpad=0.4, borderaxespad=0,
               title="RCA isoform", title_fontsize=6)
    fig.legend(handles=[Patch(color=c, label=k) for k, c in cpn_colors.items()],
               loc="lower right", bbox_to_anchor=(1 - 2 / W, 2 / H), ncol=1,
               handlelength=1, handleheight=0.7, handletextpad=0.4, borderaxespad=0,
               title="Cpn60 subunit class", title_fontsize=6)
    save(fig, out, "Fig2")

    def tips(root, cls):
        return pd.DataFrame([{"tip": l.name, "reference_absent": status.get(l.name),
                              "class": cls.get(l.name, "outgroup" if "OUTGROUP" in l.name
                                               else "")} for l in leaves(root)])
    return {"Fig2a_tips": tips(rca, rca_cls), "Fig2b_tips": tips(cpn, cpn_cls),
            "Fig2a_newick": pd.DataFrame({"newick": [open(find(build, "RCA.supported.nwk")).read().strip()]}),
            "Fig2b_newick": pd.DataFrame({"newick": [open(find(build, "CPN60.supported.nwk")).read().strip()]})}


def aln_len(path):
    n, L = 0, 0
    for line in open(path):
        if line.startswith(">"):
            if n:
                break
            n = 1
        elif n:
            L += len(line.strip())
    return L


def fig3(build, poly, out):
    W, H = 180, 92
    fig = plt.figure(figsize=(W * MM, H * MM))
    src = {}
    for i, (fam, disp, aln, col) in enumerate([("RCA", "RCA", "RCA.aln.fasta", OI["blue"]),
                                               ("CPN60", "Cpn60", "CPN60.aln.fasta",
                                                OI["orange"])]):
        sub = poly[poly.family == fam].sort_values("aln_column")
        L = aln_len(find(build, aln))
        ax = axes_mm(fig, 14, 5 + i * 45, 150, 34, H)
        ax.bar(sub.aln_column, sub.entropy_bits, width=1.0, color=col, linewidth=0)
        mean = sub.entropy_bits.mean()
        ax.axhline(mean, color="black", lw=0.6, ls=(0, (3, 2)))
        top = sub.loc[sub.entropy_bits.idxmax()]
        ax.annotate("", xy=(top.aln_column, top.entropy_bits + 0.05),
                    xytext=(top.aln_column, top.entropy_bits + 0.55),
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color="black",
                                    mutation_scale=5))
        ax.text(L * 1.01, mean, f"Mean\n{mean:.2f} bits", va="center", ha="left",
                fontsize=5.5)
        ax.set_xlim(0, L)
        ax.set_ylim(0, 3.6)
        ax.set_ylabel("Shannon entropy (bits)")
        ax.set_xlabel(f"{disp} alignment column")
        label(fig, 1, 2 + i * 45, "ab"[i], H)
        src[f"Fig3{'ab'[i]}"] = sub[["aln_column", "n_variants", "residues", "entropy_bits",
                                     "gap_fraction"]]
    save(fig, out, "Fig3")
    return src


def pangenome_partition(pav, cults):
    M = pav[cults].to_numpy(bool)
    prev = M.sum(1)
    N = len(cults)
    soft = math.ceil(0.9 * N)
    part = {"core": int((prev == N).sum()),
            "soft-core": int(((prev >= soft) & (prev < N)).sum()),
            "shell": int(((prev >= 2) & (prev < soft)).sum()),
            "cloud": int((prev == 1).sum())}
    return part, prev


def rarefaction(pav, cults, nperm=50, seed=0):
    M = pav[cults].to_numpy(bool)
    N = len(cults)
    rng = np.random.default_rng(seed)
    pan = np.zeros((nperm, N))
    core = np.zeros((nperm, N))
    idx = np.arange(N)
    for p in range(nperm):
        rng.shuffle(idx)
        cum_any = np.zeros(M.shape[0], bool)
        cum_all = np.ones(M.shape[0], bool)
        for k in range(N):
            col = M[:, idx[k]]
            cum_any |= col
            cum_all &= col
            pan[p, k] = cum_any.sum()
            core[p, k] = cum_all.sum()
    return pan, core


def fig4(gw, cults, out):
    part, prev = pangenome_partition(gw, cults)
    pan, core = rarefaction(gw, cults)
    W, H = 180, 64
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = axes_mm(fig, 16, 6, 62, 44, H)
    cats = ["core", "soft-core", "shell", "cloud"]
    vals = [part[c] for c in cats]
    tot = sum(vals)
    cols = [OI["blue"], OI["sky"], OI["orange"], OI["vermilion"]]
    ax.bar(range(4), vals, color=cols, width=0.7, linewidth=0)
    for i, v in enumerate(vals):
        ax.text(i, v + tot * 0.012, f"{v:,}\n({100 * v / tot:.0f}%)", ha="center",
                va="bottom", fontsize=5.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["Core\n(10)", "Soft-core\n(9)", "Shell\n(2–8)", "Cloud\n(1)"])
    ax.set_ylabel("Proteoform clusters (n)")
    ax.set_ylim(0, max(vals) * 1.25)
    thousands(ax)
    label(fig, 1, 2, "a", H)
    ax = axes_mm(fig, 104, 6, 72, 44, H)
    ks = np.arange(1, len(cults) + 1)
    for arr, lab, c in [(pan, "Pan-proteome", OI["blue"]), (core, "Core proteome",
                                                            OI["vermilion"])]:
        ax.fill_between(ks, arr.min(0), arr.max(0), color=c, alpha=0.18, linewidth=0)
        ax.plot(ks, arr.mean(0), "-o", color=c, ms=2.2, lw=0.8, label=lab)
    ax.set_xticks(ks)
    ax.set_xlabel("Cultivars sampled (n)")
    ax.set_ylabel("Proteoform clusters (n)")
    thousands(ax)
    ax.legend(loc="center right", handlelength=1.5)
    label(fig, 88, 2, "b", H)
    save(fig, out, "Fig4")
    acc = pd.DataFrame({"cultivars_sampled": ks, "pan_mean": pan.mean(0), "pan_min": pan.min(0),
                        "pan_max": pan.max(0), "core_mean": core.mean(0),
                        "core_min": core.min(0), "core_max": core.max(0)})
    return {"Fig4a": pd.DataFrame({"category": cats, "clusters": vals}), "Fig4b": acc}, part


def heat(ax, M, rows, cults, fs=5, extra=()):
    from matplotlib.colors import ListedColormap
    n, m = M.shape
    ax.pcolormesh(np.arange(m + 1) - 0.5, np.arange(n + 1) - 0.5, M,
                  cmap=ListedColormap([OI["absent"], OI["blue"]]), vmin=0, vmax=1,
                  edgecolors="white", linewidth=0.4)
    for j, (_, colours) in enumerate(extra):
        for yi, c in enumerate(colours):
            ax.add_patch(matplotlib.patches.Rectangle((m + 0.3 + j - 0.5, yi - 0.5), 1, 1,
                         facecolor=c, edgecolor="white", lw=0.4, clip_on=False))
    ax.set_xlim(-0.5, m - 0.5 + (len(extra) + 0.3 if extra else 0))
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(rows, fontsize=fs)
    labs = [CULT[c] for c in cults] + [lab for lab, _ in extra]
    ax.set_xticks(list(range(m)) + [m + 0.3 + j for j in range(len(extra))])
    ax.set_xticklabels(labs, rotation=90, fontsize=5.5)
    ax.tick_params(which="both", length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)


def ed1(pt, pav, nm, out):
    cults = CULT_ORDER
    comp = {r.proteotype_id: f"{r.proteotype_id} | {r.uniprot_id} | {CULT.get(r.representative_cultivar, r.representative_cultivar)}"
            for r in nm.itertuples()}

    def mat(fam, n=None):
        ids = set(pt.loc[(pt.family == fam) & pt.novel, "proteotype_id"])
        sub = pav[pav.proteotype_id.isin(ids)].copy()
        sub["prev"] = sub[cults].astype(int).sum(1)
        sub = sub.sort_values(["prev", "proteotype_id"], ascending=[False, True])
        if n:
            sub = sub.head(n)
        return sub
    a, b = mat("RCA"), mat("CPN60", 30)
    rh = 2.6
    W, H = 180, 6 + rh * len(b) + 24
    fig = plt.figure(figsize=(W * MM, H * MM))
    axa = axes_mm(fig, 50, 6, 30, rh * len(a), H)
    heat(axa, a[cults].astype(int).to_numpy(), [comp[i] for i in a.proteotype_id], cults)
    axb = axes_mm(fig, 146, 6, 30, rh * len(b), H)
    heat(axb, b[cults].astype(int).to_numpy(), [comp[i] for i in b.proteotype_id], cults)
    label(fig, 1, 2, "a", H)
    label(fig, 92, 2, "b", H)
    fig.legend(handles=[Patch(facecolor=OI["blue"], label="Present"),
                        Patch(facecolor=OI["absent"], label="Absent")], loc="lower left",
               bbox_to_anchor=(2 / W, 2 / H), ncol=2, handlelength=1, handleheight=0.7)
    save(fig, out, "ExtendedDataFig1")
    return {"EDFig1a": a.drop(columns="prev"), "EDFig1b": b.drop(columns="prev")}


def ed2(pt, pav, dom, iso, out):
    cults = CULT_ORDER
    novel = dict(zip(pt.proteotype_id, pt.novel))
    rca_col = {"RCA1\u03b2": OI["green"], "RCA2\u03b1": OI["orange"], "RCA2\u03b2": OI["sky"],
               "RCA-like": OI["purple"]}
    cpn_col = {"plastid_Cpn60alpha": OI["orange"], "plastid_Cpn60beta": OI["green"],
               "mitochondrial_CPN60": OI["purple"]}
    cls = {**{r.proteotype_id: rca_col[r.isoform] for r in iso.itertuples()},
           **{r.proteotype_id: cpn_col[r.hmm_class] for r in dom.itertuples()
              if r.hmm_class in cpn_col}}
    r = pav[pav.family == "RCA"].sort_values("proteotype_id")
    c = pav[pav.family == "CPN60"].sort_values("proteotype_id")
    blocks = [("RCA", r)] + [("Cpn60", c.iloc[i:i + 47]) for i in range(0, len(c), 47)]
    rh = 2.5
    W, H = 180, 8 + rh * 47 + 44
    fig = plt.figure(figsize=(W * MM, H * MM))
    x = 12
    for fam, sub in blocks:
        ax = axes_mm(fig, x, 8, 27, rh * len(sub), H)
        rows = [p.replace("RCA_", "").replace("CPN60_", "") for p in sub.proteotype_id]
        extra = [("Reference-absent", [OI["vermilion"] if novel[p] else "white"
                                       for p in sub.proteotype_id]),
                 ("Isoform" if fam == "RCA" else "Subunit class",
                  [cls.get(p, "white") for p in sub.proteotype_id])]
        heat(ax, sub[cults].astype(int).to_numpy(), rows, cults, extra=extra)
        ax.set_title(f"{fam} PT{sub.proteotype_id.iloc[0].split('PT')[-1]}\u2013"
                     f"PT{sub.proteotype_id.iloc[-1].split('PT')[-1]}", fontsize=6, pad=3)
        x += 43
    label(fig, 1, 2, "a", H)
    label(fig, 46, 2, "b", H)
    y0 = 2 / H
    fig.legend(handles=[Patch(facecolor=OI["blue"], label="Present"),
                        Patch(facecolor=OI["absent"], label="Absent"),
                        Patch(facecolor=OI["vermilion"], label="Reference-absent")],
               loc="lower left", bbox_to_anchor=(2 / W, y0), ncol=1, handlelength=1,
               handleheight=0.7)
    fig.legend(handles=[Patch(facecolor=v, label=k) for k, v in rca_col.items()],
               loc="lower left", bbox_to_anchor=(40 / W, y0), ncol=1, handlelength=1,
               handleheight=0.7, title="a, RCA isoform", title_fontsize=6)
    fig.legend(handles=[Patch(facecolor=cpn_col[k], label=v) for k, v in
                        [("plastid_Cpn60alpha", "Plastid Cpn60\u03b1"),
                         ("plastid_Cpn60beta", "Plastid Cpn60\u03b2"),
                         ("mitochondrial_CPN60", "Mitochondrial CPN60")]],
               loc="lower left", bbox_to_anchor=(78 / W, y0), ncol=1, handlelength=1,
               handleheight=0.7, title="b, Cpn60 subunit class", title_fontsize=6)
    save(fig, out, "ExtendedDataFig2")
    m = pav.sort_values(["family", "proteotype_id"]).copy()
    m["reference_absent"] = m.proteotype_id.map(novel)
    m["class"] = m.proteotype_id.map({**dict(zip(iso.proteotype_id, iso.isoform)),
                                      **dict(zip(dom.proteotype_id, dom.hmm_class))})
    return {"EDFig2_matrix": m}


def ed3(gw, cults, out):
    _, prev = pangenome_partition(gw, cults)
    N = len(cults)
    W, H = 180, 64
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = axes_mm(fig, 16, 6, 62, 44, H)
    counts = [int((prev == k).sum()) for k in range(1, N + 1)]
    ax.bar(range(1, N + 1), counts, color=OI["blue"], width=0.75, linewidth=0)
    ax.set_xticks(range(1, N + 1))
    ax.set_xlabel("Cultivars carrying the cluster (n)")
    ax.set_ylabel("Proteoform clusters (n)")
    thousands(ax)
    label(fig, 1, 2, "a", H)
    M = gw[cults].to_numpy(bool)
    coremask = prev == N
    cc = [int((M[:, j] & coremask).sum()) for j in range(N)]
    ac = [int((M[:, j] & ~coremask).sum()) for j in range(N)]
    order = sorted(range(N), key=lambda j: -(cc[j] + ac[j]))
    ax = axes_mm(fig, 104, 6, 72, 36, H)
    xs = np.arange(N)
    ax.bar(xs, [cc[j] for j in order], color=OI["blue"], width=0.7, label="Core", linewidth=0)
    ax.bar(xs, [ac[j] for j in order], bottom=[cc[j] for j in order], color=OI["orange"],
           width=0.7, label="Accessory", linewidth=0)
    ax.set_xticks(xs)
    ax.set_xticklabels([CULT[cults[j]] for j in order], rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel("Proteoform clusters (n)")
    thousands(ax)
    ax.legend(loc="lower right", ncol=2, handlelength=1, handleheight=0.7,
              bbox_to_anchor=(1, 1.0))
    label(fig, 88, 2, "b", H)
    save(fig, out, "ExtendedDataFig3")
    return {"EDFig3a": pd.DataFrame({"cultivars": range(1, N + 1), "clusters": counts}),
            "EDFig3b": pd.DataFrame({"cultivar": [CULT[cults[j]] for j in order],
                                     "core": [cc[j] for j in order],
                                     "accessory": [ac[j] for j in order]})}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--source-data", default=None)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    pt = pd.read_csv(find(a.inputs, "proteotypes.csv"))
    pt["novel"] = truthy(pt.novel_vs_reference)
    pav = pd.read_csv(find(a.inputs, "pav_matrix.csv"))
    poly = pd.read_csv(find(a.inputs, "polymorphic_sites.csv"))
    iso = pd.read_csv(find(a.inputs, "rca_isoform_assignment.csv"))
    dom = pd.read_csv(find(a.inputs, "domain_classification.csv"))
    nm = pd.read_csv(find(a.inputs, "proteotype_naming.csv"))
    gw = pd.read_csv(find(a.inputs, "genome_wide_pav.csv"))
    cults = [c for c in CULT_ORDER if c in gw.columns]
    src = {}
    src["Fig. 1"] = fig1(pt, pav, a.outdir)
    src["Fig. 2"] = fig2(a.inputs, pt, iso, dom, a.outdir)
    src["Fig. 3"] = fig3(a.inputs, poly, a.outdir)
    s4, part = fig4(gw, cults, a.outdir)
    src["Fig. 4"] = s4
    src["Extended Data Fig. 1"] = ed1(pt, pav, nm, a.outdir)
    src["Extended Data Fig. 2"] = ed2(pt, pav, dom, iso, a.outdir)
    src["Extended Data Fig. 3"] = ed3(gw, cults, a.outdir)
    print("[nature_figures] genome-wide partition:", part)
    if a.source_data:
        os.makedirs(a.source_data, exist_ok=True)
        for fig, sheets in src.items():
            fn = os.path.join(a.source_data, "Source_Data_" + fig.replace(". ", "").replace(" ", "_") + ".xlsx")
            with pd.ExcelWriter(fn) as xw:
                pd.DataFrame({"Source data": [f"Source data for {fig}", "Sheets: " + ", ".join(sheets)]}
                             ).to_excel(xw, sheet_name="README", index=False)
                for sheet, df in sheets.items():
                    df.to_excel(xw, sheet_name=sheet[:31], index=False)
        print("[nature_figures] source data ->", a.source_data)
    print("[nature_figures] figures ->", a.outdir)


if __name__ == "__main__":
    main()
