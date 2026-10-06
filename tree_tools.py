#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tree_tools.py - post-process and render the FastTree phylogenies of the wheat"""
from __future__ import annotations

import argparse
import math
import os
import sys
from collections import Counter

import numpy as np

def parse_newick(text):
    """Newick -> nested dict {name, length, support, children}."""
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]
    pos = 0

    def node():
        nonlocal pos
        n = {"name": None, "length": 0.0, "support": None, "children": []}
        if pos < len(text) and text[pos] == "(":
            pos += 1
            while True:
                n["children"].append(node())
                if text[pos] == ",":
                    pos += 1
                elif text[pos] == ")":
                    pos += 1
                    break
        start = pos
        while pos < len(text) and text[pos] not in ":,()":
            pos += 1
        label = text[start:pos].strip().strip("'\"")
        if n["children"]:
            try:
                n["support"] = float(label) if label else None
            except ValueError:
                n["support"] = None
        else:
            n["name"] = label
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",()":
                pos += 1
            try:
                n["length"] = float(text[start:pos])
            except ValueError:
                n["length"] = 0.0
        return n

    return node()


def to_newick(n, support_decimals=3):
    """Serialise back to Newick, keeping support labels."""
    if n["children"]:
        inner = ",".join(to_newick(c, support_decimals) for c in n["children"])
        sup = "" if n.get("support") is None else f"{n['support']:.{support_decimals}f}"
        return f"({inner}){sup}:{n['length']:.9f}"
    return f"{n['name']}:{n['length']:.9f}"


def leaves(n):
    if not n["children"]:
        return [n]
    out = []
    for c in n["children"]:
        out.extend(leaves(c))
    return out


def internal_nodes(n):
    if not n["children"]:
        return []
    out = [n]
    for c in n["children"]:
        out.extend(internal_nodes(c))
    return out


def support_summary(root):
    """Distribution of support over internal nodes that actually carry a value."""
    vals = [x["support"] for x in internal_nodes(root) if x.get("support") is not None]
    if not vals:
        return {"n": 0}
    a = np.asarray(vals, float)
    return {
        "n": int(a.size),
        "median": float(np.median(a)),
        "mean": float(a.mean()),
        "ge_0.95": int((a >= 0.95).sum()),
        "ge_0.90": int((a >= 0.90).sum()),
        "ge_0.80": int((a >= 0.80).sum()),
        "ge_0.70": int((a >= 0.70).sum()),
        "lt_0.50": int((a < 0.50).sum()),
    }


def collapse_unsupported(root, threshold):
    """Collapse internal nodes with support < threshold into their parent."""
    n_collapsed = [0]

    def walk(n):
        if not n["children"]:
            return
        for c in list(n["children"]):
            walk(c)
        merged = []
        for c in n["children"]:
            if (c["children"] and c.get("support") is not None
                    and c["support"] < threshold):
                for gc in c["children"]:
                    gc["length"] += c["length"]
                    merged.append(gc)
                n_collapsed[0] += 1
            else:
                merged.append(c)
        n["children"] = merged

    walk(root)
    return n_collapsed[0]


def _tip_distances(root):
    """{leaf_name: root-to-tip distance} and the path map needed for midpoint rooting."""
    dist, path = {}, {}

    def walk(n, d, trail):
        d2 = d + max(0.0, n["length"])
        t2 = trail + [n]
        if n["children"]:
            for c in n["children"]:
                walk(c, d2, t2)
        else:
            dist[n["name"]] = d2
            path[n["name"]] = t2

    walk(root, -max(0.0, root["length"]), [])
    return dist, path


def _reroot_on_edge(root, target, frac=0.5):
    """Re-root the tree on the branch leading to `target`, `frac` along it."""
    parents = {}

    def index(n):
        for c in n["children"]:
            parents[id(c)] = n
            index(c)

    index(root)

    chain = []
    cur = target
    while id(cur) in parents:
        chain.append(cur)
        cur = parents[id(cur)]
    old_root = cur

    blen = max(0.0, target["length"])
    new_root = {"name": None, "length": 0.0, "support": None, "children": []}
    below = target
    above_len = blen * (1.0 - frac)
    below["length"] = blen * frac

    prev = new_root
    node = parents[id(target)]
    child = target
    carry = above_len
    while True:
        node["children"] = [c for c in node["children"] if c is not child]
        nxt = parents.get(id(node))
        node_len = node["length"]
        node["length"] = carry
        if nxt is None:
            if len(node["children"]) == 1:
                only = node["children"][0]
                only["length"] += node["length"]
                prev["children"].append(only)
            else:
                prev["children"].append(node)
            break
        prev["children"].append(node)
        prev = node
        child = node
        carry = node_len
        node = nxt

    new_root["children"].insert(0, below)
    return new_root


def midpoint_root(root):
    """Root at the midpoint of the longest tip-to-tip path."""
    dist, path = _tip_distances(root)
    if len(dist) < 2:
        return root
    names = list(dist)

    def farthest_from(a):
        best, bd = None, -1.0
        pa = set(id(x) for x in path[a])
        for b in names:
            if b == a:
                continue
            mrca_d = 0.0
            for x in path[b]:
                if id(x) in pa:
                    mrca_d = max(mrca_d, dist_to_node(x))
            d = dist[a] + dist[b] - 2 * mrca_d
            if d > bd:
                best, bd = b, d
        return best, bd

    depth_cache = {}

    def dist_to_node(n):
        if id(n) in depth_cache:
            return depth_cache[id(n)]
        depth_cache[id(n)] = sum(max(0.0, x["length"]) for x in _chain_to(n))
        return depth_cache[id(n)]

    parents = {}

    def index(n):
        for c in n["children"]:
            parents[id(c)] = n
            index(c)

    index(root)

    def _chain_to(n):
        out, cur = [], n
        while id(cur) in parents:
            out.append(cur)
            cur = parents[id(cur)]
        return out

    a, _ = farthest_from(names[0])
    b, diam = farthest_from(a)
    if diam <= 0:
        return root

    half = diam / 2.0
    acc, node = 0.0, None
    cur = None
    for x in path[a][::-1]:
        pass
    cur = next(l for l in leaves(root) if l["name"] == a)
    while id(cur) in parents:
        L = max(0.0, cur["length"])
        if acc + L >= half:
            frac = (half - acc) / L if L > 0 else 0.5
            return _reroot_on_edge(root, cur, frac=1.0 - frac)
        acc += L
        cur = parents[id(cur)]
    return root


def outgroup_root(root, outgroup_names):
    """Root on the MRCA of the supplied outgroup tips."""
    want = set(outgroup_names)
    parents = {}

    def index(n):
        for c in n["children"]:
            parents[id(c)] = n
            index(c)

    index(root)
    best = None
    for n in [root] + internal_nodes(root) + leaves(root):
        names = {l["name"] for l in leaves(n)}
        if want and want <= names:
            if best is None or len(names) < len({l["name"] for l in leaves(best)}):
                best = n
    if best is None or best is root:
        sys.stderr.write("[tree_tools] outgroup not found or spans the root; "
                         "falling back to midpoint rooting.\n")
        return midpoint_root(root)
    return _reroot_on_edge(root, best, frac=0.5)


def clade_support(root, members):
    """Is `members` recovered as a clade? Return (found, support, extra, missing)."""
    want = set(members)
    best, best_names = None, None
    for n in internal_nodes(root):
        names = {l["name"] for l in leaves(n)}
        if want <= names:
            if best is None or len(names) < len(best_names):
                best, best_names = n, names
    if best is None:
        return False, None, set(), want
    return (best_names == want), best.get("support"), best_names - want, set()


def render_circular(root, out_png, novel_ids=None, prefix="", label_fontsize=None,
                    show_support_above=0.7, title_note=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    BLUE, RED = "#3B6FB6", "#C24B4B"
    novel_ids = novel_ids or set()
    tips = leaves(root)
    n = len(tips)
    if n < 3:
        return None

    gap = math.radians(18)
    span = 2 * math.pi - gap
    start = math.pi / 2 + gap / 2
    for i, lf in enumerate(tips):
        lf["angle"] = start + span * i / (n - 1)

    maxr = [0.0]

    def assign(x, r0):
        r = r0 + max(0.0, x["length"])
        x["radius"] = r
        maxr[0] = max(maxr[0], r)
        if x["children"]:
            for c in x["children"]:
                assign(c, r)
            angs = [c["angle"] for c in x["children"]]
            x["angle"] = 0.5 * (min(angs) + max(angs))

    assign(root, 0.0)
    R = maxr[0] or 1.0
    big = n > 60
    fs = label_fontsize or (3.2 if big else 6.3)

    fig, ax = plt.subplots(figsize=(14, 14) if big else (9.5, 9.5))
    ax.set_aspect("equal")
    ax.axis("off")

    def draw(x, r_par):
        a = x["angle"]
        ax.plot([r_par * math.cos(a), x["radius"] * math.cos(a)],
                [r_par * math.sin(a), x["radius"] * math.sin(a)], color="#444", lw=0.6)
        if x["children"]:
            angs = [c["angle"] for c in x["children"]]
            arc = np.linspace(min(angs), max(angs), 60)
            ax.plot(x["radius"] * np.cos(arc), x["radius"] * np.sin(arc),
                    color="#444", lw=0.6)
            s = x.get("support")
            if s is not None and s >= show_support_above and x is not root:
                ax.text(x["radius"] * math.cos(a), x["radius"] * math.sin(a),
                        f"{s:.2f}".lstrip("0"), fontsize=(2.6 if big else 4.6),
                        color="#2a6f2a", ha="center", va="center", zorder=5,
                        bbox=dict(boxstyle="round,pad=0.06", fc="white",
                                  ec="none", alpha=0.75))
            for c in x["children"]:
                draw(c, x["radius"])
        else:
            deg = math.degrees(a) % 360
            rot, ha = (deg - 180, "right") if 90 < deg < 270 else (deg, "left")
            r_lab = R * 1.06
            ax.plot([x["radius"] * math.cos(a), r_lab * math.cos(a)],
                    [x["radius"] * math.sin(a), r_lab * math.sin(a)],
                    color="#ccc", lw=0.3, zorder=0)
            is_novel = x["name"] in novel_ids
            ax.plot([x["radius"] * math.cos(a)], [x["radius"] * math.sin(a)],
                    marker="o", ms=2.2, color=(RED if is_novel else BLUE))
            ax.text((r_lab + 0.015 * R) * math.cos(a), (r_lab + 0.015 * R) * math.sin(a),
                    (x["name"] or "").replace(prefix, ""), rotation=rot,
                    rotation_mode="anchor", ha=ha, va="center", fontsize=fs,
                    color=(RED if is_novel else "#222"))

    draw(root, 0.0)
    lim = R * 1.55
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    handles = [Patch(facecolor=BLUE, label="present in Chinese Spring"),
               Patch(facecolor=RED, label="reference-absent")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.085),
               ncol=2, frameon=False, fontsize=10)
    note = title_note or (f"Node labels: SH-like local support (shown where "
                          f"≥ {show_support_above:g})")
    fig.text(0.5, 0.030, note, ha="center", fontsize=9, color="#444")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out_png


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True, help="input Newick file")
    ap.add_argument("--report", action="store_true", help="print a support summary")
    ap.add_argument("--collapse", type=float, default=None,
                    help="collapse internal nodes with support below this (e.g. 0.8)")
    ap.add_argument("--root", choices=["none", "midpoint", "outgroup"], default="none")
    ap.add_argument("--outgroup", default="", help="comma-separated tip names")
    ap.add_argument("--clade", default="",
                    help="comma-separated tips; report whether they form a clade")
    ap.add_argument("--out-nwk", default=None)
    ap.add_argument("--out-png", default=None)
    ap.add_argument("--proteotypes", default=None,
                    help="proteotypes.csv, to colour tips by reference status")
    ap.add_argument("--show-support-above", type=float, default=0.7)
    a = ap.parse_args(argv)

    root = parse_newick(open(a.tree).read())
    name = os.path.basename(a.tree)

    if a.report:
        s = support_summary(root)
        print(f"[{name}] tips={len(leaves(root))} internal nodes with support={s.get('n', 0)}")
        if s.get("n"):
            print(f"  median support {s['median']:.3f}  mean {s['mean']:.3f}")
            for k in ("ge_0.95", "ge_0.90", "ge_0.80", "ge_0.70", "lt_0.50"):
                print(f"  {k:9s}: {s[k]:4d}  ({100*s[k]/s['n']:.1f}%)")

    if a.clade:
        members = [x for x in a.clade.replace("\n", "").split(",") if x.strip()]
        mono, sup, extra, missing = clade_support(root, members)
        print(f"[{name}] clade test on {len(members)} tips: "
              f"monophyletic={mono} support={sup}")
        if extra:
            print(f"  smallest containing clade also holds: {sorted(extra)}")
        if missing:
            print(f"  tips not found in tree: {sorted(missing)}")

    if a.collapse is not None:
        k = collapse_unsupported(root, a.collapse)
        print(f"[{name}] collapsed {k} internal nodes with support < {a.collapse}")

    if a.root == "midpoint":
        root = midpoint_root(root)
        print(f"[{name}] midpoint-rooted")
    elif a.root == "outgroup":
        og = [x for x in a.outgroup.split(",") if x.strip()]
        root = outgroup_root(root, og)
        print(f"[{name}] rooted on outgroup ({len(og)} tips)")

    if a.out_nwk:
        with open(a.out_nwk, "w") as fh:
            fh.write(to_newick(root) + ";\n")
        print(f"[{name}] wrote {a.out_nwk}")

    if a.out_png:
        novel = set()
        prefix = "RCA_" if "RCA" in name.upper() else "CPN60_"
        if a.proteotypes and os.path.exists(a.proteotypes):
            import pandas as pd
            pt = pd.read_csv(a.proteotypes)
            pt["novel"] = pt["novel_vs_reference"].astype(str).str.lower().isin(["true", "1"])
            novel = set(pt.loc[pt["novel"], "proteotype_id"])
        note = "Node labels: SH-like local support"
        if a.collapse is not None:
            note += f"; nodes with support < {a.collapse:g} collapsed"
        note += {"midpoint": "; midpoint-rooted",
                 "outgroup": "; outgroup-rooted", "none": "; unrooted"}[a.root]
        render_circular(root, a.out_png, novel_ids=novel, prefix=prefix,
                        show_support_above=a.show_support_above, title_note=note)
        print(f"[{name}] wrote {a.out_png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
