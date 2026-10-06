#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""presubmission_checks.py - resolve the three open questions before submission."""
from __future__ import annotations

import argparse
import os
import re
import sys

import pandas as pd


def find(root, name):
    for c in (os.path.join(root, name),
              os.path.join(root, "audit", name),
              os.path.join(root, "domains", name),
              os.path.join(root, "transcripts", name),
              os.path.join(root, "supplementary", name)):
        if os.path.exists(c):
            return c
    return None


def read(root, name):
    p = find(root, name)
    return pd.read_csv(p) if p else None


def hdr(t):
    print("\n" + "=" * 74)
    print(f" {t}")
    print("=" * 74)


def check_pf21228(root, outdir):
    hdr("1. RCA proteotypes lacking the RCA-diagnostic domain PF21228")
    dom = read(root, "domain_classification.csv")
    pt = read(root, "proteotypes.csv")
    mem = read(root, "members.csv")
    if dom is None:
        print("  domain_classification.csv not found - run domain_profiling.py first")
        return None
    if "has_PF21228" not in dom.columns:
        print("  no PF21228 column; re-run domain_profiling.py")
        return None

    rca = dom[dom["hmm_class"] == "rubisco_activase"]
    missing = rca[~rca["has_PF21228"].astype(bool)]
    print(f"  RCA-class proteotypes            : {len(rca)}")
    print(f"  carrying PF21228                 : {int(rca['has_PF21228'].astype(bool).sum())}")
    print(f"  LACKING PF21228                  : {len(missing)}")
    if missing.empty:
        print("  Nothing to resolve.")
        return None

    rows = []
    for r in missing.itertuples():
        pid = r.proteotype_id
        info = {"proteotype_id": pid,
                "hmm_bitscore": getattr(r, "hmm_bitscore", None),
                "has_PF00004_AAA": getattr(r, "has_PF00004", None),
                "has_PF00118": getattr(r, "has_PF00118", None),
                "confidence": getattr(r, "confidence", None)}
        if pt is not None and pid in set(pt["proteotype_id"]):
            p = pt[pt["proteotype_id"] == pid].iloc[0]
            info["length_aa"] = int(p["length"])
            info["n_cultivars"] = int(p["n_cultivars"])
            info["reference_absent"] = str(p["novel_vs_reference"]).lower() in ("true", "1")
            fam = pt[pt["family"] == p["family"]]["length"]
            med = float(fam.median())
            info["family_median_length_aa"] = round(med, 1)
            info["pct_of_median_length"] = round(100 * int(p["length"]) / med, 1)
            info["verdict"] = ("TRUNCATED - a partial gene model; one clause in the text "
                               "covers it" if int(p["length"]) < 0.85 * med else
                               "FULL LENGTH - inspect before it reaches Table 5")
        if mem is not None:
            m = mem[mem["proteotype_id"] == pid]
            info["cultivars"] = ";".join(sorted(m["cultivar"].unique()))
            info["gene_loci"] = ";".join(sorted(set(m["gene"]))) if "gene" in m else ""
            info["chromosomes"] = ";".join(sorted(set(m["chromosome"]))) \
                if "chromosome" in m else ""
        cand = read(root, "candidate_table.csv")
        info["in_candidate_table_5"] = bool(cand is not None and
                                            pid in set(cand["proteotype_id"]))
        rows.append(info)

    out = pd.DataFrame(rows)
    for r in out.itertuples():
        print(f"\n  {r.proteotype_id}")
        print(f"    length            : {getattr(r,'length_aa','?')} aa "
              f"(family median {getattr(r,'family_median_length_aa','?')}, "
              f"{getattr(r,'pct_of_median_length','?')}%)")
        print(f"    other domains     : PF00004(AAA)={getattr(r,'has_PF00004_AAA','?')} "
              f"PF00118={getattr(r,'has_PF00118','?')}")
        print(f"    cultivars         : {getattr(r,'cultivars','?')}")
        print(f"    chromosome(s)     : {getattr(r,'chromosomes','?')}")
        print(f"    in Table 5        : {getattr(r,'in_candidate_table_5','?')}")
        print(f"    VERDICT           : {getattr(r,'verdict','?')}")
    out.to_csv(os.path.join(outdir, "check1_pf21228_missing.csv"), index=False)
    print(f"\n  -> check1_pf21228_missing.csv")
    return out


def check_cs_loci(root, outdir):
    hdr("2. Chinese Spring Cpn60 locus count: 31 or 32?")
    mem = read(root, "members.csv")
    if mem is None:
        print("  members.csv not found")
        return None
    if "gene" not in mem.columns:
        mem["gene"] = mem["source_id"].map(lambda s: re.sub(r"\.\d+$", "", str(s)))
    if "chromosome" not in mem.columns:
        def ch(s):
            m = re.match(r"^Traes[A-Z0-9]{2,4}(\d[ABD]|U)", str(s))
            return "un" if (not m or m.group(1) == "U") else m.group(1)
        mem["chromosome"] = mem["source_id"].map(ch)

    cs = mem[(mem["cultivar"] == "chinese_spring") & (mem["family"] == "CPN60")]
    total = cs["gene"].nunique()
    anchored = cs[cs["chromosome"] != "un"]["gene"].nunique()
    unanchored = cs[cs["chromosome"] == "un"]["gene"].nunique()
    print(f"  total distinct loci                 : {total}")
    print(f"    on named chromosomes              : {anchored}")
    print(f"    on unanchored scaffolds ('un')    : {unanchored}")
    print(f"\n  The reportable number is {total} (all loci). {anchored} counts only "
          f"chromosome-anchored loci.")
    if total != anchored:
        print(f"  The {total} vs {anchored} difference is entirely the {unanchored} "
              f"unanchored locus/loci - state which convention the text uses.")
    per_chr = (cs.groupby("chromosome")["gene"].nunique().reset_index()
               .rename(columns={"gene": "loci"}).sort_values("chromosome"))
    print("\n  per chromosome:")
    for r in per_chr.itertuples():
        print(f"    {r.chromosome:4s} {r.loci}")

    allc = (mem[mem["family"] == "CPN60"].groupby("cultivar")["gene"].nunique()
            .reset_index().rename(columns={"gene": "cpn60_loci"}))
    allr = (mem[mem["family"] == "RCA"].groupby("cultivar")["gene"].nunique()
            .reset_index().rename(columns={"gene": "rca_loci"}))
    both = allc.merge(allr, on="cultivar")
    print(f"\n  Cpn60 loci per cultivar: {both['cpn60_loci'].min()}-{both['cpn60_loci'].max()}"
          f"   RCA: {both['rca_loci'].min()}-{both['rca_loci'].max()}"
          f"{'  (INVARIANT)' if both['rca_loci'].nunique() == 1 else ''}")

    summ = pd.DataFrame([{"cultivar": "chinese_spring", "family": "CPN60",
                          "loci_total": total, "loci_anchored": anchored,
                          "loci_unanchored": unanchored}])
    summ.to_csv(os.path.join(outdir, "check2_cs_locus_count.csv"), index=False)
    per_chr.to_csv(os.path.join(outdir, "check2_cs_loci_per_chromosome.csv"), index=False)
    both.to_csv(os.path.join(outdir, "check2_loci_per_cultivar.csv"), index=False)
    print("\n  -> check2_cs_locus_count.csv, check2_cs_loci_per_chromosome.csv, "
          "check2_loci_per_cultivar.csv")
    return summ


def check_rca_sanity(root, outdir):
    hdr("3. Is the RCA prevalence discordance arithmetically forced?")
    mem = read(root, "members.csv")
    pt = read(root, "proteotypes.csv")
    prev = read(root, "prevalence_agreement.csv")
    if mem is None or pt is None:
        print("  members.csv / proteotypes.csv not found")
        return None
    if "gene" not in mem.columns:
        mem["gene"] = mem["source_id"].map(lambda s: re.sub(r"\.\d+$", "", str(s)))

    n_cult = mem["cultivar"].nunique()
    rca_m = mem[mem["family"] == "RCA"]
    per = rca_m.groupby("cultivar")["gene"].nunique()
    invariant = per.nunique() == 1
    rca_p = pt[pt["family"] == "RCA"]
    private = int((rca_p["n_cultivars"] == 1).sum())
    core = int((rca_p["n_cultivars"] == n_cult).sum())

    print(f"  cultivars                              : {n_cult}")
    print(f"  RCA loci per cultivar                  : {per.min()}"
          f"{'' if invariant else '-' + str(per.max())}"
          f"{'  (identical in every cultivar)' if invariant else ''}")
    print(f"  RCA proteotypes                        : {len(rca_p)}")
    print(f"    present in all {n_cult} cultivars (core)   : {core}")
    print(f"    cultivar-private (cloud)             : {private}")
    print()
    if invariant:
        print(f"  The nine RCA loci are present in all {n_cult} cultivars by our own count,")
        print(f"  so any orthogroup-level classification MUST call them core. Since "
              f"{private} of")
        print(f"  {len(rca_p)} RCA proteotypes are cultivar-private, at least "
              f"{private} cultivar-private")
        print(f"  proteotypes must fall in core orthogroups. The discordance is not a")
        print(f"  disagreement to explain - it is forced by the arithmetic, and its size")
        print(f"  measures allelic diversity at conserved loci.")

    ct = None
    if prev is not None and len(prev):
        p = prev.merge(pt[["proteotype_id", "family"]], on="proteotype_id", how="left")
        r = p[p["family"] == "RCA"]
        if len(r):
            ct = pd.crosstab(r["our_prevalence"], r["denovo_prevalence"])
            print("\n  RCA only - ours (rows) vs the independent annotation (columns):")
            for line in ct.to_string().split("\n"):
                print("    " + line)
            cc = int(ct.loc["cloud", "core"]) if ("cloud" in ct.index
                                                  and "core" in ct.columns) else 0
            print(f"\n  {cc} RCA proteotypes are cultivar-private by sequence yet sit in "
                  f"core orthogroups,")
            print(f"  on a family where we independently counted the loci and found them "
                  f"universal.")
            ct.to_csv(os.path.join(outdir, "check3_rca_prevalence_crosstab.csv"))
    else:
        print("\n  prevalence_agreement.csv not found - run transcript_crossref.py "
              "for the crosstab")

    out = pd.DataFrame([{"n_cultivars": n_cult,
                         "rca_loci_per_cultivar_min": int(per.min()),
                         "rca_loci_per_cultivar_max": int(per.max()),
                         "rca_loci_invariant": bool(invariant),
                         "rca_proteotypes": int(len(rca_p)),
                         "rca_core_proteotypes": core,
                         "rca_private_proteotypes": private,
                         "cloud_in_core_forced_minimum": private}])
    out.to_csv(os.path.join(outdir, "check3_rca_sanity.csv"), index=False)
    print("\n  -> check3_rca_sanity.csv"
          + (", check3_rca_prevalence_crosstab.csv" if ct is not None else ""))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", required=True, help="pan_proteome_build_* directory")
    ap.add_argument("--outdir", default=None, help="default <inputs>/checks")
    a = ap.parse_args(argv)
    root = a.inputs
    outdir = a.outdir or os.path.join(root, "checks")
    os.makedirs(outdir, exist_ok=True)

    check_pf21228(root, outdir)
    check_cs_loci(root, outdir)
    check_rca_sanity(root, outdir)

    print("\n" + "=" * 74)
    print(f" all CSVs -> {os.path.abspath(outdir)}")
    print(" upload these four for the manuscript edit:")
    for f in ("check1_pf21228_missing.csv", "check2_cs_locus_count.csv",
              "check3_rca_sanity.csv", "check3_rca_prevalence_crosstab.csv"):
        p = os.path.join(outdir, f)
        print(f"   {'x' if os.path.exists(p) else ' '} {f}")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
