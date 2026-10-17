#!/usr/bin/env python3
import argparse
import csv
import hashlib
import math
import os
import shutil
import zipfile

import pandas as pd
import yaml
from openpyxl.utils import get_column_letter

CULT = {"chinese_spring": "Chinese Spring", "arinalrfor": "ArinaLrFor", "jagger": "Jagger",
        "julius": "Julius", "lancer": "Lancer", "landmark": "Landmark", "mace": "Mace",
        "norin61": "Norin61", "stanley": "CDC Stanley", "mattis": "SY Mattis"}
FAM = {"RCA": "RCA", "CPN60": "Cpn60"}
CLASS_LABEL = {"plastid_Cpn60beta": "plastid Cpn60β", "plastid_Cpn60alpha": "plastid Cpn60α",
               "mitochondrial_CPN60": "mitochondrial CPN60"}
PROMOTER = {"reported": "Reported", "not detected": "Not detected", "not examined": "n.e."}
NA = "—"
TREES = ["RCA.aln.fasta", "CPN60.aln.fasta", "RCA.nwk", "CPN60.nwk", "RCA.supported.nwk",
         "CPN60.supported.nwk"]
CATALOGUE = ["PanProteome_catalog.annotated.faa", "WholeProteome_catalog.annotated.faa",
             "wheat_pan_proteome.db"]
SUPP_DATA = [(1, ["diamond_hits_all.csv"]), (2, ["cpn60_classes.csv"]),
             (3, ["domain_classification.csv"]), (4, ["members.csv"]), (5, ["proteotypes.csv"]),
             (6, ["pav_matrix.csv"]), (7, ["novel_divergence.csv"]),
             (8, ["rca_isoform_assignment.csv"]), (9, ["RCA_alignment.pdf"]),
             (10, ["proteotype_annotation.csv"]), (11, ["polymorphic_sites.csv"]),
             (12, ["CPN60_alignment.pdf"]), (13, ["physicochemical_features.csv"]),
             (14, ["proteotype_naming.csv"]), (15, ["proteotype_recovery.csv"]),
             (16, ["genome_wide_pav.csv"]), (17, TREES), (18, ["group_comparison.csv"])]
SUBDIRS = ["catalogue", "alignments_and_trees", "reference", "tables", "supplementary_data",
           "figures", "source_data"]


class Finder:
    def __init__(self, root, extra, exclude):
        self.root, self.extra, self.exclude = root, extra, os.path.realpath(exclude)

    def __call__(self, name):
        cands = [os.path.join(self.root, name)]
        cands += [os.path.join(self.root, d, name)
                  for d in ("audit", "domains", "transcripts", "supplementary", "figures")]
        cands += [os.path.join(d, name) for d in self.extra]
        for c in cands:
            if os.path.isfile(c):
                return c
        hits = []
        for d, _, files in os.walk(self.root):
            if os.path.realpath(d).startswith(self.exclude):
                continue
            if name in files:
                hits.append(os.path.join(d, name))
        return sorted(hits, key=lambda p: (p.count(os.sep), p))[0] if hits else None

    def table(self, name):
        p = self(name)
        if p is None:
            raise FileNotFoundError(name)
        return pd.read_csv(p, float_precision="round_trip")


def isnum(x):
    try:
        return not math.isnan(float(x))
    except (TypeError, ValueError):
        return False


def pct(v):
    s = f"{round(float(v), 1):.1f}"
    return s[:-2] if s.endswith(".0") else s


def table1(f):
    pt = f.table("proteotypes.csv")
    pt["novel"] = pt["novel_vs_reference"].astype(str).str.lower().isin(["true", "1"])
    mem = f.table("members.csv")
    sp = f.table("specificity_audit.csv").set_index("family")
    div = f.table("novel_divergence.csv")
    nq = f.table("reference_query_sets.csv").groupby("family").size()
    n_cult = int(pt["n_cultivars"].max())

    def column(fam, short_if_invariant):
        s, m = pt[pt["family"] == fam], mem[mem["family"] == fam]
        per = m.groupby("cultivar")["gene"].nunique()
        loci = (f"{per.min()}" if short_if_invariant and per.nunique() == 1
                else f"{per.min()}-{per.max()}")
        r = sp.loc[fam]
        d = div[div["family"] == fam]
        return [str(int(nq.get(fam, 0))), str(m["source_id"].nunique()), str(m["gene"].nunique()),
                loci, f"{r.identity_median:.1f}% ({r.identity_min:.1f}-{r.identity_max:.1f}%)",
                str(len(s)), str(int((s["n_cultivars"] == n_cult).sum())),
                str(int((s["n_cultivars"] == 1).sum())),
                f"{int(s['novel'].sum())} ({round(100 * float(s['novel'].mean()), 1)}%)",
                f"{round(float(d['identity_pct'].median()), 1)}%"]

    metrics = ["Curated reference queries", "Protein models recovered", "Gene loci",
               "Gene loci per cultivar", "Identity to reference query, median (range)",
               "Non-redundant proteoforms", "Core proteoforms", "Cultivar-private proteoforms",
               "Reference-absent proteoforms",
               "Median identity of reference-absent to nearest reference"]
    return [["Metric", "RCA", "Cpn60"]] + [list(r) for r in zip(metrics, column("RCA", True),
                                                                   column("CPN60", False))]


def table2(f):
    lc = f.table("locus_counts_by_cultivar.csv").sort_values("cultivar")
    cols = ["RCA_gene_loci", "RCA_protein_models", "CPN60_gene_loci", "CPN60_protein_models"]
    rows = [["Cultivar", "RCA gene loci", "RCA protein models", "Cpn60 gene loci",
             "Cpn60 protein models"]]
    rows += [[CULT.get(c, c)] + [str(int(v)) for v in vals]
             for c, vals in zip(lc["cultivar"], lc[cols].to_numpy())]
    rows.append(["Total"] + [str(int(lc[c].sum())) for c in cols])
    return rows


def table3(f):
    cs = f.table("cs_locus_table.csv")
    rows = [["Family", "Chinese Spring gene locus", "Chr", "Annotated isoforms", "Proteoform(s)"]]
    rows += [[r.family, r.gene, str(r.chromosome), str(int(r.annotated_isoforms)),
              "; ".join(str(r.proteotypes).split(";"))] for r in cs.itertuples()]
    return rows


def table4(f):
    g = f.table("granularity_partition.csv")
    rows = [["Family", "Granularity", "Clusters", "Core", "Soft-core", "Shell", "Cloud",
             "Core (%)", "Accessory (%)"]]
    for r in g.itertuples():
        counts = [r.core, r.soft_core, r.shell, r.cloud]
        rows.append([FAM.get(r.family, r.family), str(r.granularity).replace("_", " "),
                     str(int(r.clusters))]
                    + [str(int(float(x))) if isnum(x) else "N/A" for x in counts]
                    + [pct(100 * float(r.core) / float(r.clusters)) if isnum(r.core) else "N/A",
                       pct(r.accessory_pct) if isnum(r.accessory_pct) else "N/A"])
    return rows


def table5(f, n_cult):
    s = f.table("rca_isoform_summary.csv")
    rows = [["Established isoform", "Chinese Spring locus", "Proteoforms (n)", "Full-length (n)",
             "Reference-absent (n)", "Principal full-length proteoform (cultivars)",
             f"Cultivars with a full-length form (of {n_cult})",
             "Promoter heat-responsive element"]]
    for r in s.itertuples():
        el = str(r.promoter_heat_responsive_element)
        rows.append([r.established_isoform, r.chinese_spring_locus, str(int(r.proteoforms)),
                     str(int(r.full_length)), str(int(r.reference_absent)),
                     f"{r.principal_full_length_proteoform} ({int(r.principal_n_cultivars)})",
                     str(int(r.cultivars_with_full_length_form)), PROMOTER.get(el, el)])
    return rows


def table6(f):
    p = f.table("physicochemical_features.csv")
    c = f.table("cpn60_classes.csv")
    spec = [("Length (aa)", "length", lambda v: f"{v:g}"),
            ("Molecular weight (kDa)", "mw", lambda v: f"{v / 1000:.1f}"),
            ("Isoelectric point (pI)", "pI", lambda v: f"{v:.2f}"),
            ("GRAVY", "gravy", lambda v: f"{v:.3f}"),
            ("Aliphatic index", "aliphatic_index", lambda v: str(round(v, 2))),
            ("Instability index", "instability_index", lambda v: f"{v:.2f}")]
    rows = [["Descriptor", "RCA", "Cpn60"]]
    for label, col, fmt in spec:
        rows.append([label] + [fmt(float(p.loc[p["family"] == fam, col].median()))
                               for fam in ("RCA", "CPN60")])
    m = p.merge(c[["proteotype_id", "cpn60_class"]], on="proteotype_id")
    for k, g in m[m["family"] == "CPN60"].groupby("cpn60_class"):
        rows.append([f"Aliphatic index {str(k).replace('_', ' ')} (n={len(g)})", "-",
                     f"{g['aliphatic_index'].median():.1f}"])
    return rows


def table7(f):
    cand = f.table("candidate_table.csv")
    iso = f.table("rca_isoform_assignment.csv")
    dc = f.table("domain_classification.csv")
    label = dict(zip(iso["proteotype_id"], iso["established_isoform"]))
    for pid, cls in zip(dc["proteotype_id"], dc["hmm_class"]):
        if str(pid).startswith("CPN60_"):
            label.setdefault(pid, CLASS_LABEL.get(cls, cls))
    cand["_tier"] = (~cand["candidate_class"].astype(str).str.startswith("reference-absent")).astype(int)
    cand = cand.sort_values(["_tier", "proteotype_id"])
    rows = [["Composite identifier", "Priority tier", "Established isoform or class", "Chr",
             "Length (aa)", "Cultivars (n)", "Nearest CS proteoform", "Identity (%)",
             "Residue differences"]]
    for r in cand.itertuples():
        tier = str(r.candidate_class)
        rows.append([r.composite_id, tier[:1].upper() + tier[1:], label.get(r.proteotype_id, ""),
                     str(r.chromosome), str(int(r.length)), str(int(r.n_cultivars)),
                     str(r.nearest_reference_proteotype) if isinstance(r.nearest_reference_proteotype, str)
                     and r.nearest_reference_proteotype else NA,
                     str(float(r.identity_pct)) if isnum(r.identity_pct) else NA,
                     str(int(float(r.residue_differences))) if isnum(r.residue_differences) else NA])
    return rows


def write_rows(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(rows)


def write_xlsx(path, sheets):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets:
            df.to_excel(xw, sheet_name=name, index=False)
            ws = xw.book[name]
            ws.freeze_panes = "A2"
            for j, col in enumerate(df.columns, 1):
                w = max(10, len(str(col)) + 2)
                if len(df):
                    w = max(w, int(df[col].astype(str).str.len().quantile(0.9)) + 2)
                ws.column_dimensions[get_column_letter(j)].width = min(60, w)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build(root, outdir, cfg=None, log=None, figures=None, source_data=None):
    say = log.info if log else print
    cfg = cfg or {}
    inp = os.path.abspath((cfg.get("project") or {}).get("inputs_dir", "inputs"))
    extra = [inp, os.path.join(inp, "references"),
             os.path.abspath(os.path.join(root, "..", "..", "inputs")),
             os.path.abspath(os.path.join(root, "..", "..", "inputs", "references"))]
    os.makedirs(outdir, exist_ok=True)
    f = Finder(root, extra, outdir)
    for d in SUBDIRS:
        shutil.rmtree(os.path.join(outdir, d), ignore_errors=True)
        os.makedirs(os.path.join(outdir, d))
    for s in ("main", "extended_data"):
        os.makedirs(os.path.join(outdir, "figures", s))
    missing = []

    def copy(name, sub, dst_name=None):
        p = f(name)
        if p is None:
            missing.append(name)
            return None
        dst = os.path.join(outdir, sub, dst_name or name)
        shutil.copyfile(p, dst)
        return dst

    for name in CATALOGUE:
        copy(name, "catalogue")
    for name in TREES:
        copy(name, "alignments_and_trees")
    copy("reference_query_sets.csv", "reference")

    n_cult = 10
    try:
        n_cult = int(f.table("proteotypes.csv")["n_cultivars"].max())
    except FileNotFoundError:
        pass
    makers = [table1, table2, table3, table4, lambda x: table5(x, n_cult), table6, table7]
    for i, maker in enumerate(makers, 1):
        try:
            write_rows(os.path.join(outdir, "tables", f"Supplementary_Table_{i}.csv"), maker(f))
        except (FileNotFoundError, KeyError) as e:
            missing.append(f"Supplementary Table {i} ({e})")

    sd = os.path.join(outdir, "supplementary_data")
    for n, names in SUPP_DATA:
        base = os.path.join(sd, f"Supplementary_Data_{n}")
        if n == 17:
            paths = [f(x) for x in names]
            if None in paths:
                missing.append(f"Supplementary Data {n}")
                continue
            with zipfile.ZipFile(base + ".zip", "w", zipfile.ZIP_DEFLATED) as z:
                for name, p in zip(names, paths):
                    z.write(p, name)
            continue
        if names[0].endswith(".pdf"):
            copy(names[0], "supplementary_data", f"Supplementary_Data_{n}.pdf")
            continue
        stem = os.path.splitext(names[0])[0]
        if copy(names[0], "tables", f"Supplementary_Data_{n}_{stem}.csv") is None:
            continue
        sheets = [("Data", f.table(names[0]))]
        if n == 1 and copy("reference_query_sets.csv", "tables",
                           "Supplementary_Data_1_reference_queries.csv"):
            sheets.append(("Reference queries", f.table("reference_query_sets.csv")))
        write_xlsx(base + ".xlsx", sheets)

    figdir = figures or os.path.join(root, "figures", "nature")
    if os.path.isdir(figdir):
        for name in sorted(os.listdir(figdir)):
            if name.endswith((".pdf", ".png")) and name.startswith(("Fig", "ExtendedDataFig")):
                sub = "extended_data" if name.startswith("Extended") else "main"
                shutil.copyfile(os.path.join(figdir, name), os.path.join(outdir, "figures", sub, name))
    else:
        missing.append(figdir)
    srcdir = source_data or os.path.join(root, "source_data")
    if os.path.isdir(srcdir):
        for name in sorted(os.listdir(srcdir)):
            if name.endswith(".xlsx"):
                shutil.copyfile(os.path.join(srcdir, name), os.path.join(outdir, "source_data", name))
    else:
        missing.append(srcdir)

    rows = []
    for d, _, files in os.walk(outdir):
        for name in files:
            p = os.path.join(d, name)
            rel = os.path.relpath(p, outdir)
            if rel != "MANIFEST.tsv":
                rows.append((rel, os.path.getsize(p), md5(p)))
    with open(os.path.join(outdir, "MANIFEST.tsv"), "w") as fh:
        fh.write("path\tbytes\tmd5\n")
        for r in sorted(rows):
            fh.write("\t".join(map(str, r)) + "\n")
    say(f"[build_tables] {len(rows)} files -> {outdir}")
    for m in missing:
        say(f"[build_tables] WARNING missing: {m}")
    return outdir


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--config", default=None)
    ap.add_argument("--figures", default=None)
    ap.add_argument("--source-data", default=None)
    a = ap.parse_args(argv)
    cfg = {}
    if a.config and os.path.exists(a.config):
        with open(a.config) as fh:
            cfg = yaml.safe_load(fh) or {}
    build(a.inputs, a.outdir or os.path.join(a.inputs, "final"), cfg,
          figures=a.figures, source_data=a.source_data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
