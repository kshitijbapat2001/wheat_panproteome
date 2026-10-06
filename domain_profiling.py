#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""domain_profiling.py - resolve family and subunit identity by HMM profiling and"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request

import pandas as pd

UA = {"User-Agent": "wheat-panproteome/3.1"}
INTERPRO_HMM = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}/?annotation=hmm"
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

PFAM = {
    "PF00118": "Cpn60_TCP1 (group I/II chaperonin)",
    "PF00004": "AAA (AAA+ ATPase core)",
    "PF21228": "RuBisCO_activase_AAA_helical (RCA-diagnostic)",
}

SEEDS = {
    "plastid_Cpn60alpha": '(taxonomy_id:33090) AND (reviewed:true) AND '
                          '(protein_name:"chaperonin 60 subunit alpha" OR '
                          'protein_name:"RuBisCO large subunit-binding protein subunit alpha")',
    "plastid_Cpn60beta": '(taxonomy_id:33090) AND (reviewed:true) AND '
                         '(protein_name:"chaperonin 60 subunit beta" OR '
                         'protein_name:"RuBisCO large subunit-binding protein subunit beta")',
    "mitochondrial_CPN60": '(reviewed:true) AND ('
                           'protein_name:"chaperonin CPN60, mitochondrial" OR '
                           'protein_name:"60 kDa chaperonin, mitochondrial" OR '
                           'protein_name:"60 kDa heat shock protein, mitochondrial" OR '
                           'gene:HSPD1 OR gene:HSP60)',
    "rubisco_activase": '(taxonomy_id:33090) AND (reviewed:true) AND '
                        '(protein_name:"ribulose bisphosphate carboxylase/oxygenase activase")',
    "cbbX_negative_control": '(reviewed:true) AND (gene:cbbX)',
}

MIN_SEEDS = 3
DEFAULT_MARGIN = 20.0


def sh(cmd, log=print, check=True):
    log("    $ " + " ".join(str(c) for c in cmd))
    p = subprocess.run([str(c) for c in cmd], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    if check and p.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed ({p.returncode})\n"
                           + "\n".join((p.stdout or "").splitlines()[-15:]))
    return p


def fetch(url, dest, log=print, binary=True):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        log(f"    cache {os.path.basename(dest)}")
        return dest
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as fh:
        shutil.copyfileobj(r, fh)
    log(f"    got   {os.path.basename(dest)} ({os.path.getsize(dest)/1024:.0f} kB)")
    return dest


def get_pfam_hmms(cache, log=print):
    """Download the Pfam HMMs and concatenate them into one pressed database."""
    combined = os.path.join(cache, "pfam_selected.hmm")
    if not (os.path.exists(combined) and os.path.getsize(combined) > 0):
        parts = []
        for acc in PFAM:
            raw = os.path.join(cache, f"{acc}.hmm.gz")
            try:
                fetch(INTERPRO_HMM.format(acc=acc), raw, log)
            except Exception as e:
                log(f"    WARNING could not fetch {acc}: {e}")
                continue
            out = raw[:-3]
            try:
                with gzip.open(raw, "rb") as gi, open(out, "wb") as go:
                    shutil.copyfileobj(gi, go)
            except OSError:
                shutil.copyfile(raw, out)
            parts.append(out)
        if not parts:
            return None
        with open(combined, "wb") as fo:
            for p in parts:
                with open(p, "rb") as fi:
                    shutil.copyfileobj(fi, fo)
                fo.write(b"\n")
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        if not os.path.exists(combined + ext):
            sh(["hmmpress", "-f", combined], log)
            break
    return combined


def get_seed_sets(cache, log=print):
    """Reviewed sequences per class from UniProt, one FASTA each."""
    out = {}
    for cls, query in SEEDS.items():
        fa = os.path.join(cache, f"seed_{cls}.faa")
        url = (UNIPROT_STREAM + "?format=fasta&size=500&query="
               + urllib.parse.quote(query))
        try:
            fetch(url, fa, log)
        except Exception as e:
            log(f"    WARNING seed fetch failed for {cls}: {e}")
            continue
        n = sum(1 for l in open(fa) if l.startswith(">"))
        log(f"    {cls:26s} {n} reviewed sequence(s)")
        if n >= MIN_SEEDS:
            out[cls] = fa
        else:
            log(f"      too few (<{MIN_SEEDS}); class profile skipped")
    return out


def build_class_hmms(seeds, work, log=print):
    """MAFFT-align each seed set and hmmbuild a profile; press them together."""
    hmms = []
    for cls, fa in seeds.items():
        aln = os.path.join(work, f"{cls}.aln")
        hmm = os.path.join(work, f"{cls}.hmm")
        if not os.path.exists(hmm):
            with open(aln, "w") as fo:
                p = subprocess.run(["mafft", "--auto", "--quiet", fa],
                                   stdout=fo, stderr=subprocess.PIPE, text=True)
            if p.returncode != 0:
                log(f"    WARNING mafft failed for {cls}"); continue
            sh(["hmmbuild", "--informat", "afa", "-n", cls, hmm, aln], log)
        hmms.append(hmm)
    if not hmms:
        return None
    combined = os.path.join(work, "class_profiles.hmm")
    with open(combined, "wb") as fo:
        for h in hmms:
            with open(h, "rb") as fi:
                shutil.copyfileobj(fi, fo)
            fo.write(b"\n")
    sh(["hmmpress", "-f", combined], log)
    return combined


def parse_domtblout(path):
    rows = []
    if not os.path.exists(path):
        return pd.DataFrame(columns=["query", "profile", "bitscore", "evalue"])
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 14:
                continue
            rows.append({"query": f[3], "profile": f[0], "accession": f[1],
                         "evalue": float(f[6]), "bitscore": float(f[7])})
    return pd.DataFrame(rows)


def run(inputs, outdir, cache=None, margin=DEFAULT_MARGIN, threads=8, log=print):
    os.makedirs(outdir, exist_ok=True)
    cache = cache or os.path.join(outdir, "cache")
    os.makedirs(cache, exist_ok=True)

    fa = None
    for c in ("PanProteome_catalog.faa", "PanProteome_catalog.annotated.faa"):
        p = os.path.join(inputs, c)
        if os.path.exists(p):
            fa = p; break
    if fa is None:
        raise SystemExit(f"no proteotype FASTA in {inputs}")
    log(f"[domains] proteotypes: {os.path.basename(fa)}")

    for tool in ("hmmscan", "hmmbuild", "hmmpress", "mafft", "diamond"):
        if not shutil.which(tool):
            raise SystemExit(f"{tool} not on PATH - activate the conda environment")

    log("[domains] fetching Pfam profiles")
    pfam = get_pfam_hmms(cache, log)
    dom = pd.DataFrame()
    if pfam:
        tbl = os.path.join(outdir, "pfam_domains.domtblout")
        if not os.path.exists(tbl):
            sh(["hmmscan", "--cpu", str(threads), "--domtblout", tbl,
                "-E", "1e-5", "-o", os.devnull, pfam, fa], log)
        dom = parse_domtblout(tbl)
        if len(dom):
            best = dom.sort_values("bitscore", ascending=False)
            wide = (best.groupby(["query", "profile"])["bitscore"].max()
                    .unstack(fill_value=0.0).reset_index())
            wide.to_csv(os.path.join(outdir, "pfam_domain_matrix.csv"), index=False)
            log(f"[domains] Pfam hits for {best['query'].nunique()} proteotypes")

    log("[domains] retrieving curated seed sets")
    seeds = get_seed_sets(cache, log)
    cls_hmm = build_class_hmms(seeds, cache, log) if seeds else None
    calls = pd.DataFrame()
    if cls_hmm:
        tbl2 = os.path.join(outdir, "class_profiles.domtblout")
        if not os.path.exists(tbl2):
            sh(["hmmscan", "--cpu", str(threads), "--domtblout", tbl2,
                "-E", "1e-3", "-o", os.devnull, cls_hmm, fa], log)
        d2 = parse_domtblout(tbl2)
        if len(d2):
            d2 = d2.sort_values("bitscore", ascending=False)
            rows = []
            for q, g in d2.groupby("query", sort=False):
                g = g.drop_duplicates("profile")
                top = g.iloc[0]
                second = g.iloc[1] if len(g) > 1 else None
                rows.append({
                    "proteotype_id": q,
                    "hmm_class": top["profile"],
                    "hmm_bitscore": round(float(top["bitscore"]), 1),
                    "hmm_evalue": top["evalue"],
                    "runner_up": (second["profile"] if second is not None else ""),
                    "bitscore_margin": round(float(top["bitscore"]) -
                                             float(second["bitscore"]), 1)
                    if second is not None else float("inf"),
                })
            calls = pd.DataFrame(rows)

    rbh = pd.DataFrame()
    if seeds:
        allseed = os.path.join(cache, "curated_all.faa")
        with open(allseed, "w") as fo:
            for cls, s in seeds.items():
                for line in open(s):
                    fo.write(line.replace(">", f">{cls}|", 1) if line.startswith(">")
                             else line)
        db1, db2 = os.path.join(cache, "seed_db"), os.path.join(cache, "pt_db")
        sh(["diamond", "makedb", "--in", allseed, "--db", db1, "--quiet"], log)
        sh(["diamond", "makedb", "--in", fa, "--db", db2, "--quiet"], log)
        f1 = os.path.join(outdir, "pt_vs_seed.tsv")
        f2 = os.path.join(outdir, "seed_vs_pt.tsv")
        cols = ["qseqid", "sseqid", "pident", "evalue", "bitscore"]
        for q, db, out in ((fa, db1, f1), (allseed, db2, f2)):
            if not os.path.exists(out):
                sh(["diamond", "blastp", "--query", q, "--db", db + ".dmnd",
                    "--out", out, "--outfmt", "6", *cols, "--max-target-seqs", "1",
                    "--max-hsps", "1", "--evalue", "1e-5", "--more-sensitive",
                    "--threads", str(threads), "--quiet"], log)
        a = pd.read_csv(f1, sep="\t", names=cols)
        b = pd.read_csv(f2, sep="\t", names=cols)
        fwd = dict(zip(a.qseqid, a.sseqid))
        rev = dict(zip(b.qseqid, b.sseqid))
        rows = []
        for pt, seed in fwd.items():
            rows.append({"proteotype_id": pt, "rbh_partner": seed,
                         "rbh_class": seed.split("|")[0],
                         "reciprocal": rev.get(seed) == pt})
        rbh = pd.DataFrame(rows)

    out = calls if len(calls) else pd.DataFrame(columns=["proteotype_id"])
    if len(rbh):
        out = out.merge(rbh, on="proteotype_id", how="outer")
    if len(dom):
        for acc, name in PFAM.items():
            hit = set(dom.loc[dom["accession"].str.startswith(acc), "query"]) | \
                  set(dom.loc[dom["profile"].str.contains(acc, na=False), "query"])
            out[f"has_{acc}"] = out["proteotype_id"].isin(hit)

    def verdict(r):
        """Confidence in the CLASS call."""
        h, rb = r.get("hmm_class"), r.get("rbh_class")
        if pd.isna(h) or not isinstance(h, str):
            return "unassigned"
        if h == "cbbX_negative_control":
            return "OUT-OF-FAMILY (CbbX-like) - inspect"
        same_class = isinstance(rb, str) and rb == h
        wide = float(r.get("bitscore_margin", 0) or 0) >= margin
        if same_class and wide:
            return "high confidence"
        if same_class:
            return "medium (class agrees, margin below threshold)"
        if pd.isna(rb) or not isinstance(rb, str):
            return "medium (profile only, no best hit)"
        return "low (profile and best-hit class disagree)"

    if len(out):
        out["confidence"] = out.apply(verdict, axis=1)
        out = out.sort_values(["hmm_class", "proteotype_id"])
        out.to_csv(os.path.join(outdir, "domain_classification.csv"), index=False)

        log("\n[domains] class assignment by HMM profile:")
        for k, v in out["hmm_class"].value_counts().items():
            log(f"    {k:28s} {v:4d}")
        log("[domains] confidence:")
        for k, v in out["confidence"].value_counts().items():
            log(f"    {k:46s} {v:4d}")
        n_cbbx = int((out["hmm_class"] == "cbbX_negative_control").sum())
        log(f"\n[domains] negative control: {n_cbbx} proteotype(s) scored best against "
            f"CbbX. Zero is the expected result in a green-plant proteome; anything "
            f"above zero needs manual inspection.")
        if "has_PF21228" in out:
            rca = out[out["hmm_class"] == "rubisco_activase"]
            if len(rca):
                log(f"[domains] of {len(rca)} RCA-class proteotypes, "
                    f"{int(rca['has_PF21228'].sum())} carry PF21228 "
                    f"(the RCA-diagnostic helical domain).")
    log(f"\n[domains] -> {os.path.abspath(outdir)}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", required=True, help="pan_proteome_build_* directory")
    ap.add_argument("--outdir", default=None, help="default <inputs>/domains")
    ap.add_argument("--cache", default=None, help="persistent download cache")
    ap.add_argument("--margin", type=float, default=DEFAULT_MARGIN)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    run(a.inputs, a.outdir or os.path.join(a.inputs, "domains"),
        a.cache, a.margin, a.threads)
    return 0


if __name__ == "__main__":
    sys.exit(main())
