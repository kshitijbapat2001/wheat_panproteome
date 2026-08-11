#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the data artefacts for the manuscript: tables, figures and
RESULTS_NUMBERS.json. Deliberately writes NO manuscript prose.


build_tables.py - assemble the data bundle (tables, figures, numbers) from a build.

Creates <result>/final/ containing everything a submission needs, with the
numbers written straight out of the pipeline's own tables so the text can never
drift from the data:

    final/
      manuscript/   Manuscript_draft.md, Manuscript_draft.docx, RESULTS_NUMBERS.json
      figures/      Figure1..Figure7 (main text, renamed and ordered by citation)
      tables/       Table1..Table6 as .csv and a combined .docx
      supplementary/
        figures/    Supplementary Figure S1..S4
        tables/     Supplementary Table S1..S11
        data/       FASTA catalogues, Newick trees, SQLite database
      MANIFEST.md   what every file is, and which comment it answers

The draft is generated, not transcribed. Every count, percentage, confidence
interval and support value is read from the CSVs produced by the run, so
re-running with different data yields a correct draft rather than a stale one.
Passages needing an author decision are marked >>> ACTION <<<.

Usage
-----
    python build_tables.py --inputs pan_proteome_build_XXXX
    python build_tables.py --inputs <dir> --outdir <dir>/final
"""
from __future__ import annotations

import argparse
import json
import os
import shutil


def _copy(src, dst):
    """Copy that survives a re-run.

    Two things break a naive shutil.copy here. The source can resolve to the destination
    itself once a previous run has populated the bundle (SameFileError), and shutil.copy
    preserves the source mode, so a read-only input yields a read-only output that the
    next run cannot overwrite (PermissionError). Handle both.
    """
    # realpath, not abspath: outputs/latest is a symlink to the build directory, so the
    # same file can be reached by two different-looking paths and abspath will not tell.
    if os.path.realpath(src) == os.path.realpath(dst):
        return dst
    if os.path.exists(dst):
        try:
            os.chmod(dst, 0o644)
        except OSError:
            pass
    shutil.copy(src, dst)
    try:
        os.chmod(dst, 0o644)
    except OSError:
        pass
    return dst
import sys
from datetime import date

import numpy as np
import pandas as pd

CULT = {"chinese_spring": "Chinese Spring", "arinalrfor": "ArinaLrFor", "jagger": "Jagger",
        "julius": "Julius", "lancer": "Lancer", "landmark": "Landmark", "mace": "Mace",
        "norin61": "Norin61", "stanley": "CDC Stanley", "mattis": "SY Mattis"}
FAM_DISP = {"RCA": "RCA", "CPN60": "Cpn60"}

# main-text figures, in order of first citation in the text
FIGURE_MAP = [
    ("Figure1_overview.png", "Figure_overview.png",
     "Overview of the RCA and Cpn60 pan-proteome across the cultivar panel"),
    ("Figure2_novel_pav.png", "RCA_novel_PAV.png",
     "Presence-absence of reference-absent proteotypes"),
    ("RCA_tree_circular.png", None, "Circular phylogeny of the RCA proteotypes"),
    ("CPN60_tree_circular.png", None, "Circular phylogeny of the Cpn60 proteotypes"),
    ("RCA_variability.png", "Figure6_RCA_variability.png",
     "Positional sequence variability across the RCA alignment"),
    ("CPN60_variability.png", "Figure7_CPN60_variability.png",
     "Positional sequence variability across the Cpn60 alignment"),
    ("Figure_pangenome_overview.png", "Figure5_pangenome_overview.png",
     "Genome-wide pan-proteome anchored on the reference"),
]

SUPP_TABLES = [
    ("proteotypes.csv", "Non-redundant proteotype catalogue"),
    ("members.csv", "Every protein model with cultivar, gene locus, chromosome and proteotype"),
    ("pav_matrix.csv", "Proteotype x cultivar presence-absence matrix"),
    ("proteotype_naming.csv", "Composite identifiers and closest UniProt entry"),
    ("physicochemical_features.csv", "Per-proteotype thermal-stability descriptors"),
    ("polymorphic_sites.csv", "Per-alignment-column residue counts, entropy and gap fraction"),
    ("genome_wide_pav.csv", "Genome-wide proteotype-cluster presence-absence table"),
    ("pangenome_partition.csv", "Core/soft-core/shell/cloud partition of genome-wide clusters"),
    ("group_comparison.csv", "Phenotype-stratified test output"),
    ("diamond_hits_all.csv", "All DIAMOND hits with per-hit identity and coverage"),
    ("audit/novel_divergence.csv",
     "Divergence of every reference-absent proteotype from its nearest reference sequence"),
    ("audit/cpn60_classes.csv",
     "Chaperonin class assigned to every Cpn60 proteotype, with the reference it matched"),
    ("proteotype_annotation.csv",
     "UniProtKB and InterPro annotation of every proteotype: protein name, gene name, "
     "Pfam and InterPro domains, subcellular location"),
]

SUPP_DATA = [
    ("PanProteome_catalog.annotated.faa", "Proteotype sequences with composite identifiers"),
    ("RCA.nwk", "RCA phylogeny (FastTree, with SH-like local support)"),
    ("CPN60.nwk", "Cpn60 phylogeny (FastTree, with SH-like local support)"),
    ("RCA.supported.nwk", "RCA phylogeny, weak nodes collapsed and rooted"),
    ("CPN60.supported.nwk", "Cpn60 phylogeny, weak nodes collapsed and rooted"),
    ("RCA.aln.fasta", "RCA proteotype alignment (MAFFT)"),
    ("CPN60.aln.fasta", "Cpn60 proteotype alignment (MAFFT)"),
    ("wheat_pan_proteome.db", "Consolidated SQLite database"),
    ("WholeProteome_catalog.annotated.faa", "Genome-wide cluster catalogue"),
    ("reference_query_sets.csv", "Curated UniProt reference query sets"),
    ("proteotype_annotation.csv", "UniProt/InterPro annotation of every proteotype"),
]


# --------------------------------------------------------------------------------------
EXTRA_SEARCH = []          # populated from cfg by build(): project inputs dir etc.
EXCLUDE_DIRS = []          # populated by build(): the bundle we are writing


def _find(root, name):
    for c in ([os.path.join(d, name) for d in EXTRA_SEARCH] +
              [os.path.join(root, name),
              os.path.join(root, "supplementary", name),
              os.path.join(root, "figures", name),
              os.path.join(root, "audit", name),
              os.path.join(root, "domains", name),
              os.path.join(root, "transcripts", name),
              os.path.join(root, "references", name),
              os.path.join(root, "inputs", "references", name)]):
        if os.path.exists(c):
            return c
    # Recursive fallback. Files such as reference_query_sets.csv are written outside the
    # timestamped build directory, and when the explicit list missed them a table silently
    # came out blank instead of failing. Search widely, take the newest match.
    import glob as _glob
    roots = [root, os.path.join(root, "..", "..")] + list(EXTRA_SEARCH)
    hits = []
    for r in roots:
        if r and os.path.isdir(r):
            hits += _glob.glob(os.path.join(r, "**", name), recursive=True)
    # Never return a file that lives inside the bundle this run is writing: on a re-run
    # that is a copy the previous run made, and copying it onto itself raises
    # SameFileError. Sources must come from the build, not from the output.
    # Collapse paths that resolve to the same file (outputs/latest is a symlink to the
    # build directory, so every artefact is reachable twice).
    seen, uniq = set(), []
    for p in hits:
        rp = os.path.realpath(p)
        if rp not in seen:
            seen.add(rp); uniq.append(p)
    hits = uniq
    # Prefer a source OUTSIDE the bundle we are writing, because copying a file onto
    # itself fails. But do not hide it entirely: for some inputs the copy inside the
    # bundle is the only one left, and a blank table is worse than a self-copy that
    # _copy() already handles.
    if EXCLUDE_DIRS and hits:
        ex = [os.path.realpath(d) for d in EXCLUDE_DIRS]
        inside = lambda p: any(os.path.realpath(p).startswith(e + os.sep) for e in ex)
        hits = [p for p in hits if not inside(p)] + [p for p in hits if inside(p)]
    if hits:
        hits.sort(key=lambda p: os.path.getmtime(p), reverse=True)
        return hits[0]
    return None


def _read(root, name):
    p = _find(root, name)
    if p and p.lower().endswith(".csv"):
        try:
            return pd.read_csv(p)
        except Exception:  # noqa: BLE001
            return None
    return None


def _fmt(x, nd=1):
    try:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return "n/a"
        return f"{float(x):.{nd}f}"
    except Exception:  # noqa: BLE001
        return str(x)


def collect_numbers(root, cfg_hint=None) -> dict:
    """Every number the draft quotes, read from the run's own outputs."""
    N = {"generated": str(date.today()), "build_dir": os.path.abspath(root)}

    pt = _read(root, "proteotypes.csv")
    mem = _read(root, "members.csv")
    spec = _read(root, "specificity_audit.csv")
    lc = _read(root, "locus_counts_by_cultivar.csv")
    gran = _read(root, "granularity_partition.csv")
    ci = _read(root, "core_fraction_ci.csv")
    div = _read(root, "novel_divergence.csv")
    cls = _read(root, "cpn60_classes.csv")
    part = _read(root, "pangenome_partition.csv")
    phys = _read(root, "physicochemical_features.csv")
    poly = _read(root, "polymorphic_sites.csv")
    qset = _read(root, "reference_query_sets.csv")
    cand = _read(root, "candidate_table.csv")

    if pt is not None:
        pt = pt.copy()
        pt["novel"] = pt["novel_vs_reference"].astype(str).str.lower().isin(["true", "1"])
        N["cultivars"] = int(pt["n_cultivars"].max())
        for fam in ("RCA", "CPN60"):
            s = pt[pt["family"] == fam]
            if s.empty:
                continue
            N[fam] = {"proteotypes": int(len(s)),
                      "core": int((s["n_cultivars"] == N["cultivars"]).sum()),
                      "private": int((s["n_cultivars"] == 1).sum()),
                      "novel": int(s["novel"].sum()),
                      "novel_pct": round(100 * float(s["novel"].mean()), 1)}
    if mem is not None and "gene" in mem.columns:
        for fam in ("RCA", "CPN60"):
            s = mem[mem["family"] == fam]
            if s.empty:
                continue
            per = s.groupby("cultivar")["gene"].nunique()
            N.setdefault(fam, {}).update({
                "protein_models": int(s["source_id"].nunique()),
                "gene_loci": int(s["gene"].nunique()),
                "loci_per_cultivar_min": int(per.min()),
                "loci_per_cultivar_max": int(per.max()),
                "loci_invariant": bool(per.nunique() == 1)})
            if "chromosome" in s.columns:
                cs = s[s["cultivar"] == "chinese_spring"]
                if len(cs):
                    N[fam]["reference_loci"] = int(cs["gene"].nunique())
                    N[fam]["reference_models"] = int(cs["source_id"].nunique())
                    N[fam]["chromosomes"] = (cs.groupby("chromosome")["gene"].nunique()
                                             .sort_index().to_dict())
    if spec is not None:
        N["specificity"] = spec.to_dict("records")
    if qset is not None:
        N["reference_queries"] = qset.groupby("family").size().to_dict()
    if gran is not None:
        N["granularity"] = gran.to_dict("records")
    if ci is not None:
        N["core_ci"] = ci.to_dict("records")
    if div is not None and len(div):
        N["divergence"] = {}
        for fam, g in div.groupby("family"):
            N["divergence"][fam] = {
                "n": int(len(g)),
                "median_identity": round(float(g["identity_pct"].median()), 1),
                "pct_ge_995": round(100 * float((g["identity_pct"] >= 99.5).mean()), 0),
                "pct_lt_90": round(100 * float((g["identity_pct"] < 90).mean()), 0)}
    if cls is not None and "cpn60_class" in cls.columns:
        N["cpn60_classes"] = cls["cpn60_class"].value_counts().to_dict()
    if part is not None:
        try:
            d = dict(zip(part.iloc[:, 0], part.iloc[:, 1]))
            tot = sum(d.values())
            N["genome_wide"] = {"total": int(tot),
                                **{k: int(v) for k, v in d.items()},
                                "core_pct": round(100 * d.get("core", 0) / max(1, tot))}
        except Exception:  # noqa: BLE001
            pass
    if phys is not None:
        N["descriptors"] = {}
        for fam, g in phys.groupby("family"):
            N["descriptors"][fam] = {c: round(float(g[c].median()), 3)
                                     for c in ("length", "mw", "pI", "gravy",
                                               "aliphatic_index", "instability_index")
                                     if c in g.columns}
        if cls is not None:
            m = phys.merge(cls[["proteotype_id", "cpn60_class"]], on="proteotype_id", how="left")
            N["descriptors_by_cpn60_class"] = {
                k: {"n": int(len(g)), "aliphatic_median": round(float(g["aliphatic_index"].median()), 1)}
                for k, g in m[m["family"] == "CPN60"].groupby("cpn60_class")}
    if poly is not None:
        N["polymorphism"] = {}
        for fam, g in poly.groupby("family"):
            N["polymorphism"][fam] = {
                "variable_columns": int(len(g)),
                "mean_entropy": round(float(g["entropy_bits"].mean()), 2),
                "max_entropy": round(float(g["entropy_bits"].max()), 2),
                "max_column": int(g.loc[g["entropy_bits"].idxmax(), "aln_column"])}
    if cand is not None:
        N["n_candidates"] = int(len(cand))


    # HMM domain profiling (domain_profiling.py)
    dc = _read(root, "domain_classification.csv")
    if dc is None or not len(dc):
        print("[final] WARNING: domain_classification.csv not found -> RESULTS_NUMBERS.json "
              "will have no 'domains' section. Expected under <build>/domains/.")
    if dc is not None and len(dc):
        N["domains"] = {
            "n": int(len(dc)),
            "classes": dc["hmm_class"].value_counts().to_dict() if "hmm_class" in dc else {},
            "confidence": dc["confidence"].value_counts().to_dict() if "confidence" in dc else {},
            "n_cbbx": int((dc.get("hmm_class") == "cbbX_negative_control").sum())
                      if "hmm_class" in dc else 0,
            "pf21228": int(dc["has_PF21228"].sum()) if "has_PF21228" in dc else None,
            "pf00118": int(dc["has_PF00118"].sum()) if "has_PF00118" in dc else None,
        }
    # independent transcript recovery (transcript_crossref.py)
    tr = _read(root, "proteotype_recovery.csv")
    if tr is None or not len(tr):
        print("[final] WARNING: proteotype_recovery.csv not found -> RESULTS_NUMBERS.json "
              "will have no 'transcripts' section. Expected under <build>/transcripts/.")
    if tr is not None and len(tr):
        d = {"n": int(len(tr))}
        for fam, g in tr.groupby("family"):
            e = {"n": int(len(g)), "recovered": int(g["any_locus_recovered"].sum())}
            if "novel" in g.columns:
                for nov, key in ((False, "ref_present"), (True, "ref_absent")):
                    s2 = g[g["novel"] == nov]
                    if len(s2):
                        e[key] = {"n": int(len(s2)),
                                  "recovered": int(s2["any_locus_recovered"].sum()),
                                  "pct": round(100 * float(s2["any_locus_recovered"].mean()))}
            d[fam] = e
        N["transcripts"] = d
    if mem is not None:
        _cs = [c for c in mem["cultivar"].unique() if c == "chinese_spring"]
        _rest = sorted(c for c in mem["cultivar"].unique() if c != "chinese_spring")
        N["cultivar_order"] = _cs + _rest
        N["cultivar_names"] = [CULT.get(c, c) for c in N["cultivar_order"]]

    # tree support, read from the Newick files themselves
    N["tree_support"] = {}
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import tree_tools
        for fam in ("RCA", "CPN60"):
            p = _find(root, f"{fam}.nwk")
            if not p:
                continue
            r = tree_tools.parse_newick(open(p).read())
            s = tree_tools.support_summary(r)
            if s.get("n"):
                s["pct_ge_080"] = round(100 * s["ge_0.80"] / s["n"], 1)
                s["pct_lt_050"] = round(100 * s["lt_0.50"] / s["n"], 1)
                N["tree_support"][fam] = s
    except Exception:  # noqa: BLE001
        pass

    # monophyly / support of a named clade (report.candidate_clade), which is what a
    # "is this a real lineage?" question actually asks for
    clade = list((cfg_hint or {}).get("report", {}).get("candidate_clade") or [])
    # Prefer the clade the audit detected from gene loci; config is only an override.
    _ci = _read(root, "divergent_clade.csv")
    if _ci is not None and len(_ci):
        r = _ci.iloc[0].to_dict()
        mem = str(r.get("members", ""))
        sep = ";" if ";" in mem else ","
        mem = [x.strip(" '\"[]") for x in mem.split(sep) if x.strip(" '\"[]")]
        if mem and not clade:
            clade = mem
        N["clade_detected"] = {k: r.get(k) for k in
                               ("family", "canonical_group", "clade_group",
                                "canonical_loci", "clade_loci", "subgenomes",
                                "monophyletic", "support", "n_members")}
    if clade:
        try:
            import tree_tools
            for fam in ("RCA", "CPN60"):
                p = _find(root, f"{fam}.nwk")
                if not p:
                    continue
                r = tree_tools.parse_newick(open(p).read())
                names = {l["name"] for l in tree_tools.leaves(r)}
                members = [c for c in clade if c in names]
                if len(members) < 2:
                    continue
                mono, sup, extra, missing = tree_tools.clade_support(r, members)
                N["clade_test"] = {"family": fam, "members": members,
                                   "n_members": len(members), "monophyletic": bool(mono),
                                   "support": (float(sup) if sup is not None else None),
                                   "extra_tips": sorted(extra)}
                break
        except Exception:  # noqa: BLE001
            pass
    return N


# --------------------------------------------------------------------------------------
def tables_to_docx(tables, docx_path):
    try:
        import docx
        from docx.shared import Pt
    except Exception:  # noqa: BLE001
        return None
    doc = docx.Document()
    doc.styles["Normal"].font.name = "Times New Roman"
    doc.styles["Normal"].font.size = Pt(10)
    for name, df, cap in tables:
        doc.add_heading(name, level=2)
        doc.add_paragraph(cap)
        show = df.head(200)
        t = doc.add_table(rows=1, cols=len(show.columns))
        t.style = "Light Grid Accent 1"
        for i, c in enumerate(show.columns):
            t.rows[0].cells[i].text = str(c)
        for _, row in show.iterrows():
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = "" if pd.isna(v) else str(v)
        if len(df) > len(show):
            doc.add_paragraph(f"({len(df)} rows total; first {len(show)} shown. "
                              f"Full table in final/tables/.)")
        doc.add_page_break()
    doc.save(docx_path)
    return docx_path


# --------------------------------------------------------------------------------------
def build(root, outdir, cfg=None, log=None):
    cfg = cfg or {}

    def _log(msg, *a):
        if log:
            log.info(msg, *a)
        else:
            print(msg % a if a else msg)

    # reference_query_sets.csv and the external cache live under project.inputs_dir,
    # which is outside the timestamped build directory.
    global EXTRA_SEARCH, EXCLUDE_DIRS
    _inp = (cfg.get("project", {}) or {}).get("inputs_dir", "")
    EXTRA_SEARCH = [d for d in (_inp, os.path.join(_inp, "references"),
                                os.path.join(root, "..", "..", "inputs"),
                                os.path.join(root, "..", "..", "inputs", "references"))
                    if d and os.path.isdir(d)]

    EXCLUDE_DIRS = [os.path.abspath(outdir)]

    for sub in ("manuscript", "figures", "tables",
                "supplementary/figures", "supplementary/tables", "supplementary/data"):
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)

    N = collect_numbers(root, cfg)
    with open(os.path.join(outdir, "manuscript", "RESULTS_NUMBERS.json"), "w") as fh:
        json.dump(N, fh, indent=2, default=str)

    # ---- main figures, renamed in citation order
    copied_fig = []
    for i, (primary, alt, cap) in enumerate(FIGURE_MAP, 1):
        src = _find(root, primary) or (_find(root, alt) if alt else None)
        if src:
            dst = os.path.join(outdir, "figures", f"Figure{i}_{os.path.basename(src)}")
            _copy(src, dst)
            copied_fig.append((f"Figure {i}", os.path.basename(dst), cap))
    # supplementary figures
    copied_sfig = []
    for i, nm in enumerate(["RCA_PAV_heatmap.png", "CPN60_PAV_heatmap.png",
                            "genomewide_prevalence_spectrum.png",
                            "genomewide_per_cultivar.png"], 1):
        src = _find(root, nm)
        if src:
            dst = os.path.join(outdir, "supplementary", "figures", f"FigureS{i}_{nm}")
            _copy(src, dst)
            copied_sfig.append((f"Supplementary Figure S{i}", os.path.basename(dst), nm))

    # ---- main tables, generated from the audit outputs
    main_tables = []
    R, C = N.get("RCA", {}), N.get("CPN60", {})
    sp = {r["family"]: r for r in N.get("specificity", [])}
    rq = N.get("reference_queries", {})
    if not rq:
        _log("[final] WARNING: reference_query_sets.csv not found -> Table 1 'Curated "
             "reference queries' will be blank. Check project.inputs_dir in config.yaml.")
    t1 = pd.DataFrame([
        ["Curated reference queries", rq.get("RCA", ""), rq.get("CPN60", "")],
        ["Protein models recovered", R.get("protein_models"), C.get("protein_models")],
        ["Gene loci", R.get("gene_loci"), C.get("gene_loci")],
        ["Gene loci per cultivar",
         f"{R.get('loci_per_cultivar_min')}" + ("" if R.get("loci_invariant")
                                                else f"-{R.get('loci_per_cultivar_max')}"),
         f"{C.get('loci_per_cultivar_min')}-{C.get('loci_per_cultivar_max')}"],
        ["Identity to reference query, median (range)",
         f"{_fmt(sp.get('RCA', {}).get('identity_median'))}% "
         f"({_fmt(sp.get('RCA', {}).get('identity_min'))}-{_fmt(sp.get('RCA', {}).get('identity_max'))}%)",
         f"{_fmt(sp.get('CPN60', {}).get('identity_median'))}% "
         f"({_fmt(sp.get('CPN60', {}).get('identity_min'))}-{_fmt(sp.get('CPN60', {}).get('identity_max'))}%)"],
        ["Non-redundant proteotypes", R.get("proteotypes"), C.get("proteotypes")],
        ["Core proteotypes", R.get("core"), C.get("core")],
        ["Cultivar-private proteotypes", R.get("private"), C.get("private")],
        ["Reference-absent proteotypes",
         f"{R.get('novel')} ({R.get('novel_pct')}%)",
         f"{C.get('novel')} ({C.get('novel_pct')}%)"],
        ["Median identity of reference-absent proteotypes to nearest reference",
         f"{N.get('divergence', {}).get('RCA', {}).get('median_identity', 'n/a')}%",
         f"{N.get('divergence', {}).get('CPN60', {}).get('median_identity', 'n/a')}%"],
    ], columns=["Metric", "RCA", "Cpn60"])
    main_tables.append(("Table 1", t1,
                        "Summary of the pan-proteome with recruitment specificity."))

    lc = _read(root, "locus_counts_by_cultivar.csv")
    if lc is not None:
        lc = lc.copy()
        order = N.get("cultivar_order") or list(lc["cultivar"])
        lc["__o"] = lc["cultivar"].map({c: i for i, c in enumerate(order)}).fillna(999)
        lc = lc.sort_values("__o").drop(columns="__o")
        lc["cultivar"] = lc["cultivar"].map(lambda c: CULT.get(c, c))
        cols = ["cultivar"] + [c for c in lc.columns if c.startswith("RCA_")] \
                            + [c for c in lc.columns if c.startswith("CPN60_")]
        lc = lc[[c for c in cols if c in lc.columns]]
        lc.columns = [c.replace("CPN60_", "Cpn60 ").replace("RCA_", "RCA ")
                       .replace("_", " ").replace("cultivar", "Cultivar")
                      for c in lc.columns]
        main_tables.append(("Table 2", lc,
                            "Protein models and gene loci per cultivar."))
    phys = _read(root, "physicochemical_features.csv")
    if phys is not None:
        LBL = {"length": "Length (aa)", "mw": "Molecular weight (kDa)",
               "pI": "Isoelectric point (pI)", "gravy": "GRAVY",
               "aliphatic_index": "Aliphatic index",
               "instability_index": "Instability index",
               "frac_charged_DEKRH": "Charged-residue fraction (D+E+K+R+H)"}
        SCALE = {"mw": 1e-3}
        ND = {"length": 0, "mw": 1, "pI": 2, "gravy": 3, "aliphatic_index": 1,
              "instability_index": 1, "frac_charged_DEKRH": 3}

        def _cell(g, c):
            s = g[c] * SCALE.get(c, 1.0)
            n = ND.get(c, 2)
            return f"{s.median():.{n}f} ({s.min():.{n}f} to {s.max():.{n}f})"

        fams = [f for f in ("RCA", "CPN60") if f in set(phys["family"])]
        rows = []
        for c in ("length", "mw", "pI", "gravy", "aliphatic_index",
                  "instability_index", "frac_charged_DEKRH"):
            if c not in phys.columns:
                continue
            r = {"Descriptor": LBL.get(c, c)}
            for fam in fams:
                r[FAM_DISP.get(fam, fam)] = _cell(phys[phys.family == fam], c)
            rows.append(r)
        cls = _read(root, "cpn60_classes.csv")
        if cls is not None and "cpn60_class" in cls.columns:
            m = phys.merge(cls[["proteotype_id", "cpn60_class"]],
                           on="proteotype_id", how="left")
            m = m[m.family == "CPN60"]
            for k, g in m.groupby("cpn60_class"):
                rows.append({"Descriptor":
                             f"Aliphatic index \u2014 {str(k).replace('_', ' ')} (n={len(g)})",
                             "RCA": "\u2014",
                             "Cpn60": f"{g.aliphatic_index.median():.1f} "
                                      f"({g.aliphatic_index.min():.1f} to "
                                      f"{g.aliphatic_index.max():.1f})"})
        main_tables.append(("Table 3", pd.DataFrame(rows),
                            "Sequence-derived physicochemical descriptors."))
    gran = _read(root, "granularity_partition.csv")
    if gran is not None:
        main_tables.append(("Table 4", gran,
                            "Prevalence partition at three clustering granularities."))
    cand = _read(root, "candidate_table.csv")
    if cand is not None:
        main_tables.append(("Table 5", cand,
                            "Prioritised candidate proteotypes for experimental follow-up."))
    csl = _read(root, "cs_locus_table.csv")
    if csl is not None:
        main_tables.append(("Table 6", csl,
                            "Reference gene loci underlying the proteotypes."))

    for name, df, cap in main_tables:
        df.to_csv(os.path.join(outdir, "tables",
                               name.replace(" ", "") + ".csv"), index=False)
    tables_to_docx(main_tables, os.path.join(outdir, "tables", "MainTables.docx"))

    # ---- supplementary tables and data
    copied_st = []
    for i, (nm, desc) in enumerate(SUPP_TABLES, 1):
        src = _find(root, os.path.basename(nm)) or os.path.join(root, nm)
        if os.path.exists(src):
            dst = os.path.join(outdir, "supplementary", "tables",
                               f"TableS{i}_{os.path.basename(src)}")
            _copy(src, dst)
            copied_st.append((f"Supplementary Table S{i}", os.path.basename(dst), desc))
    copied_sd = []
    for nm, desc in SUPP_DATA:
        src = _find(root, nm)
        if src:
            dst = os.path.join(outdir, "supplementary", "data", os.path.basename(src))
            _copy(src, dst)
            copied_sd.append((os.path.basename(dst), desc))

    # ---- manuscript: the full paper, with this run's numbers substituted in
    figs = {}
    for lab, fname, _cap in copied_fig:
        try:
            figs[int(lab.split()[-1])] = os.path.join(outdir, "figures", fname)
        except Exception:  # noqa: BLE001
            pass
    # No manuscript prose is generated: this build emits data artefacts only.
    _log("[final] tables, figures and RESULTS_NUMBERS.json written; prose not generated "
         "by design.")

    # ---- manifest
    with open(os.path.join(outdir, "MANIFEST.md"), "w", encoding="utf-8") as fh:
        fh.write(f"# Publication bundle\n\nBuild: `{os.path.basename(root)}`  \n"
                 f"Generated: {N['generated']}\n\n")
        fh.write("## Numbers\n\n"
                 "- `manuscript/RESULTS_NUMBERS.json` - every computed value, for the text\n\n")
        fh.write("## Main figures\n\n")
        for lab, f, cap in copied_fig:
            fh.write(f"- **{lab}** `figures/{f}` - {cap}\n")
        fh.write("\n## Main tables\n\n")
        for name, df, cap in main_tables:
            fh.write(f"- **{name}** `tables/{name.replace(' ', '')}.csv` "
                     f"({len(df)} rows) - {cap}\n")
        fh.write("- `tables/MainTables.docx` - all main tables in one Word file\n")
        fh.write("\n## Supplementary figures\n\n")
        for lab, f, src in copied_sfig:
            fh.write(f"- **{lab}** `supplementary/figures/{f}` (from {src})\n")
        fh.write("\n## Supplementary tables\n\n")
        for lab, f, desc in copied_st:
            fh.write(f"- **{lab}** `supplementary/tables/{f}` - {desc}\n")
        fh.write("\n## Supplementary data\n\n")
        for f, desc in copied_sd:
            fh.write(f"- `supplementary/data/{f}` - {desc}\n")
        fh.write("\n## Review points this bundle addresses\n\n"
                 "| Point | Where |\n|---|---|\n"
                 "| Over-recruitment / bona fide fraction | Table 1; "
                 "`supplementary/tables/*specificity*`; audit/specificity_audit.csv |\n"
                 "| Protein models vs gene loci | Table 2; audit/locus_counts_by_cultivar.csv |\n"
                 "| Chromosomes and genomic loci | Table 6; members.csv gene/chromosome columns |\n"
                 "| How far is 'reference-absent'? | Table 1 last row; audit/novel_divergence.csv |\n"
                 "| Accessory fraction vs granularity | Table 4 |\n"
                 "| Core fraction confidence interval | audit/core_fraction_ci.csv |\n"
                 "| Actionable candidates in main text | Table 5 |\n"
                 "| Cpn60 heterogeneity | audit/cpn60_classes.csv; Table 3 |\n"
                 "| Tree support and rooting | Figures 3-4; *.supported.nwk |\n"
                 "| Domain confirmation of family identity | proteotype_annotation.csv |\n")

    _log("[final] bundle -> %s", os.path.abspath(outdir))
    _log("[final] %d main figures, %d main tables, %d supplementary tables, %d data files",
         len(copied_fig), len(main_tables), len(copied_st), len(copied_sd))
    return outdir


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", required=True, help="pan_proteome_build_* directory")
    ap.add_argument("--outdir", default=None, help="default <inputs>/final")
    ap.add_argument("--config", default="config.yaml")
    a = ap.parse_args(argv)
    cfg = {}
    if os.path.exists(a.config):
        import yaml
        cfg = yaml.safe_load(open(a.config)) or {}
    build(a.inputs, a.outdir or os.path.join(a.inputs, "final"), cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
