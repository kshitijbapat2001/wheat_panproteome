#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
panproteome_audit.py - post-hoc audit and gene-level resolution for the wheat
RCA/Cpn60 pan-proteome pipeline.

Motivation
----------
`hpc_pipeline.py` reports *protein models* (one Ensembl record per annotated
transcript) and defines proteotypes by exact sequence identity. Both choices are
defensible, but on their own they produce three numbers that are easy to
misread and that reviewers challenged:

  1. "191 RCA copies" looks like massive over-recruitment against ~6-9 expected
     loci. It is not: it is 191 transcript isoforms of 90 genes, nine per
     cultivar. Only a gene-level resolution shows this.
  2. "52% of proteotypes are absent from the reference" sounds like presence /
     absence variation. Because a proteotype is an exact-identity cluster, a
     single residue - including an annotation error - creates one. The distance
     to the nearest reference sequence has to be measured.
  3. The Cpn60 family was searched with only TWO UniProt reference sequences at
     40% identity, which recruits the whole group I chaperonin superfamily
     including mitochondrial CPN60. Family-level statistics are meaningless
     until that is stratified.

This module computes exactly those missing quantities from the pipeline's own
outputs, writes them as supplementary tables, and returns them for the report.
It is deliberately standalone: it can be run against any completed build
directory without re-running DIAMOND.

Usage
-----
    python panproteome_audit.py --inputs pan_proteome_build_XXXX \
                                --outdir  pan_proteome_build_XXXX/audit

Outputs
-------
    gene_locus_table.csv        every protein model -> gene, chromosome, subgenome
    locus_counts_by_cultivar.csv    protein models and gene loci per cultivar
    cs_locus_table.csv          Chinese Spring loci (manuscript Table 6)
    specificity_audit.csv       identity of each model to its best reference query
    cpn60_classes.csv           per-proteotype chaperonin class assignment
    novel_divergence.csv        reference-absent proteotype -> nearest CS sequence
    granularity_partition.csv   core/shell/cloud at three clustering granularities
    core_fraction_ci.csv        bootstrap CIs on core fractions
    candidate_table.csv         prioritised candidates (manuscript Table 5)

Dependencies: pandas, numpy; biopython only for --divergence (pairwise aligner).
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys

import numpy as np
import pandas as pd

FAMILIES = ("RCA", "CPN60")

# Ensembl Plants wheat gene ids look like TraesCS4A02G177500.1 (Chinese Spring)
# or TraesJUL4A03G02120760.2 (cultivar assemblies). The 2-4 character block after
# "Traes" is the assembly code; the chromosome follows it.
_GENE_RE = re.compile(r"^(Traes[A-Z0-9]{2,4}(?:\d[ABDU]|U)\d*G\d+)\.\d+$")
_CHR_RE = re.compile(r"^Traes[A-Z0-9]{2,4}(\d[ABD]|U)")


# --------------------------------------------------------------------------------------
# 1. gene-locus resolution
# --------------------------------------------------------------------------------------
def gene_of(source_id: str) -> str:
    """Protein model id -> parent gene id (strip the transcript suffix)."""
    m = _GENE_RE.match(str(source_id))
    return m.group(1) if m else re.sub(r"\.\d+$", "", str(source_id))


def chrom_of(source_id: str) -> str:
    """Protein model id -> chromosome, or 'un' for unanchored scaffolds."""
    m = _CHR_RE.match(str(source_id))
    if not m:
        return "un"
    return "un" if m.group(1) == "U" else m.group(1)


def resolve_gene_loci(members: pd.DataFrame) -> pd.DataFrame:
    """Add gene / chromosome / subgenome / homoeologous-group columns to members."""
    out = members.copy()
    out["gene"] = out["source_id"].map(gene_of)
    out["chromosome"] = out["source_id"].map(chrom_of)
    out["subgenome"] = out["chromosome"].str[-1].where(out["chromosome"] != "un", "un")
    out["homoeologous_group"] = out["chromosome"].str[0].where(out["chromosome"] != "un", "un")
    return out


def locus_counts(members_g: pd.DataFrame) -> pd.DataFrame:
    """Protein models and distinct gene loci per cultivar per family."""
    g = (members_g.groupby(["family", "cultivar"])
         .agg(protein_models=("source_id", "nunique"), gene_loci=("gene", "nunique"))
         .reset_index())
    piv = g.pivot(index="cultivar", columns="family",
                  values=["protein_models", "gene_loci"])
    piv.columns = [f"{fam}_{metric}" for metric, fam in piv.columns]
    return piv.reset_index()


def cs_locus_table(members_g: pd.DataFrame, reference="chinese_spring") -> pd.DataFrame:
    """Reference loci with isoform counts and the proteotypes they define."""
    cs = members_g[members_g["cultivar"] == reference]
    t = (cs.groupby(["family", "gene", "chromosome"])
         .agg(annotated_isoforms=("source_id", "nunique"),
              proteotypes=("proteotype_id", lambda s: ";".join(sorted(set(s)))))
         .reset_index()
         .sort_values(["family", "chromosome", "gene"]))
    return t


# --------------------------------------------------------------------------------------
# 2. recruitment specificity
# --------------------------------------------------------------------------------------
def parse_ref_acc(qseqid: str) -> str:
    parts = str(qseqid).split("|")
    return parts[1] if len(parts) >= 3 else str(qseqid)


def specificity_audit(hits: pd.DataFrame, members_g: pd.DataFrame,
                      qid2family: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Best identity of every retained protein model to its family's reference set.

    Returns (per-model table, per-family summary). A tight, high-identity
    distribution means the search recruited the intended family; a broad or
    multimodal one means it recruited a superfamily and the results must be
    stratified before they are interpreted.
    """
    h = hits.copy()
    h["ref_acc"] = h["qseqid"].map(parse_ref_acc)
    fam_of_model = dict(zip(members_g["source_id"], members_g["family"]))
    h = h[h["sseqid"].isin(fam_of_model)]
    h["family"] = h["sseqid"].map(fam_of_model)

    best = (h.sort_values("bitscore", ascending=False)
            .groupby(["family", "sseqid"], as_index=False).head(1)
            [["family", "cultivar", "sseqid", "ref_acc", "pident", "qcovhsp",
              "scovhsp", "evalue", "bitscore"]]
            .rename(columns={"sseqid": "source_id", "pident": "identity_to_reference"}))

    rows = []
    for fam, g in best.groupby("family"):
        a = g["identity_to_reference"].to_numpy(float)
        bands = {}
        for lo, hi in [(40, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 95), (95, 101)]:
            bands[f"n_{lo}_{hi}"] = int(((a >= lo) & (a < hi)).sum())
        rows.append({"family": fam, "n_models": int(a.size),
                     "n_reference_queries": int(h.loc[h.family == fam, "ref_acc"].nunique()),
                     "identity_min": round(float(a.min()), 1),
                     "identity_median": round(float(np.median(a)), 1),
                     "identity_mean": round(float(a.mean()), 1),
                     "identity_max": round(float(a.max()), 1), **bands})
    return best, pd.DataFrame(rows)


def _class_from_description(desc: str) -> str | None:
    """Map a UniProt FASTA description to a chaperonin class."""
    d = (desc or "").lower()
    if not d:
        return None
    if "mitochondri" in d:
        return "mitochondrial_groupI"
    if "alpha" in d or "subunit alpha" in d:
        return "plastid_Cpn60alpha"
    if "beta" in d or "subunit-binding protein beta" in d:
        return "plastid_Cpn60beta"
    if "chloroplast" in d or "plastid" in d or "rubisco" in d:
        return "plastid_chaperonin_unspecified"
    if "chaperonin" in d or "cpn60" in d or "hsp60" in d or "groel" in d:
        return "chaperonin_unassigned_compartment"
    return None


def assign_cpn60_classes(best: pd.DataFrame, members_g: pd.DataFrame,
                         beta_ref="P08823", class2_ref="A0A2I4JD18",
                         beta_cut=85.0, class2_cut=77.0, intermediate_cut=55.0,
                         hits: pd.DataFrame | None = None,
                         query_table: str | None = None,
                         specific_margin: float = 15.0) -> pd.DataFrame:
    """Provisional chaperonin class per proteotype, by identity to the two queries.

    Cut-offs are empirical, taken from the modes of the identity distribution.
    This is a stopgap: definitive assignment needs HMM domain profiling and
    reciprocal best hits against characterised alpha/beta subunits. It is
    reported as provisional in the manuscript for exactly that reason.
    """
    if hits is None:
        raise ValueError("assign_cpn60_classes needs the full hits table")

    # Preferred route: classify by what the best-matching reference actually IS.
    # reference_query_sets.csv carries the UniProt description of every query, so the
    # class comes from annotation rather than from an identity cut-off calibrated on a
    # two-sequence query set.
    ann = {}
    if query_table and os.path.exists(query_table):
        try:
            q = pd.read_csv(query_table)
            for r in q[q["family"] == "CPN60"].itertuples():
                c = _class_from_description(getattr(r, "description", ""))
                if c:
                    ann[str(r.accession)] = c
        except Exception:  # noqa: BLE001
            ann = {}
    if ann:
        hb = hits.copy()
        hb["ref_acc"] = hb["qseqid"].map(parse_ref_acc)
        hb = hb[hb["ref_acc"].isin(ann)]
        if len(hb):
            hb["cls"] = hb["ref_acc"].map(ann)
            # A reference annotated only as "Chaperonin 60-like protein" carries no
            # compartment information. If a reference with a SPECIFIC annotation
            # (mitochondrial / alpha / beta) matches within `margin` identity points of
            # the best hit, prefer it: a slightly lower identity to a well-annotated
            # reference is more informative than a marginally better hit to a vague one.
            GENERIC = {"chaperonin_unassigned_compartment",
                       "plastid_chaperonin_unspecified"}
            margin = float(specific_margin)
            hb = hb.sort_values("pident", ascending=False)
            rows = []
            for sid, g in hb.groupby("sseqid", sort=False):
                bestp = float(g["pident"].iloc[0])
                spec = g[(~g["cls"].isin(GENERIC)) & (g["pident"] >= bestp - margin)]
                pick = spec.iloc[0] if len(spec) else g.iloc[0]
                rows.append({"sseqid": sid, "ref_acc": pick["ref_acc"],
                             "pident": float(pick["pident"]), "cls": pick["cls"]})
            top = pd.DataFrame(rows)
            mm = members_g[members_g["family"] == "CPN60"].merge(
                top, left_on="source_id", right_on="sseqid", how="left")
            per = (mm.groupby("proteotype_id")
                   .agg(cpn60_class=("cls", lambda s: s.dropna().mode().iloc[0]
                                     if len(s.dropna().mode()) else "unclassified"),
                        best_reference=("ref_acc", lambda s: s.dropna().iloc[0]
                                        if len(s.dropna()) else ""),
                        identity_to_best_reference=("pident", "max"))
                   .reset_index())
            print(f"[audit] Cpn60 classes assigned from reference annotation "
                  f"({len(ann)} annotated queries): {len(per)} proteotypes")
            n_un = int((per["cpn60_class"] == "unclassified").sum())
            if n_un:
                sys.stderr.write(f"[audit] {n_un} Cpn60 proteotype(s) unclassified by "
                                 f"annotation; falling back to identity bands for those\n")
            return per

    # Fallback: identity to the two canonical references, as before.
    h = hits.copy()
    h["ref_acc"] = h["qseqid"].map(parse_ref_acc)
    h = h[h["ref_acc"].isin([beta_ref, class2_ref])]
    piv = (h.sort_values("bitscore", ascending=False)
           .groupby(["sseqid", "ref_acc"], as_index=False).head(1)
           .pivot_table(index="sseqid", columns="ref_acc", values="pident"))
    for col in (beta_ref, class2_ref):
        if col not in piv.columns:
            piv[col] = np.nan
    piv["best"] = piv[[beta_ref, class2_ref]].max(axis=1)

    def cls(r):
        if pd.isna(r["best"]):
            return "unclassified"
        if r[beta_ref] >= beta_cut:
            return "plastid_Cpn60beta_like"
        if r[class2_ref] >= class2_cut:
            return "plastid_chaperonin_classII"
        if r["best"] >= intermediate_cut:
            return "intermediate_divergent"
        return "mitochondrial_or_related_groupI"

    piv["cpn60_class"] = piv.apply(cls, axis=1)
    mm = members_g[members_g["family"] == "CPN60"].merge(
        piv.reset_index(), left_on="source_id", right_on="sseqid", how="left")
    per_pt = (mm.groupby("proteotype_id")
              .agg(cpn60_class=("cpn60_class", lambda s: s.mode().iloc[0] if len(s.mode()) else np.nan),
                   identity_to_beta=(beta_ref, "max"),
                   identity_to_classII=(class2_ref, "max"))
              .reset_index())
    return per_pt


# --------------------------------------------------------------------------------------
# 3. how far is "reference-absent"?
# --------------------------------------------------------------------------------------
def novel_divergence(proteotypes: pd.DataFrame, sequences: dict) -> pd.DataFrame:
    """Pairwise identity of every reference-absent proteotype to the nearest
    reference-present proteotype of the same family.

    Answers the objection that a proteotype is reference-absent as soon as it
    differs by one residue. Global Needleman-Wunsch with BLOSUM62.
    """
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices
    except Exception as e:  # pragma: no cover
        raise RuntimeError("Biopython is required for novel_divergence()") from e

    al = Align.PairwiseAligner()
    al.mode = "global"
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")

    def ident(a, b):
        aln = al.align(a, b)[0]
        A, B = aln[0], aln[1]
        match = sum(1 for x, y in zip(A, B) if x == y and x != "-")
        ncol = sum(1 for x, y in zip(A, B) if not (x == "-" and y == "-"))
        diff = sum(1 for x, y in zip(A, B) if x != y)
        return (100.0 * match / ncol if ncol else float("nan")), diff

    pt = proteotypes.copy()
    pt["novel"] = pt["novel_vs_reference"].astype(str).str.lower().isin(["true", "1"])
    rows = []
    for fam, g in pt.groupby("family"):
        ref_ids = [p for p in g.loc[~g["novel"], "proteotype_id"] if p in sequences]
        for r in g[g["novel"]].itertuples():
            if r.proteotype_id not in sequences or not ref_ids:
                continue
            best = None
            for c in ref_ids:
                pid, diff = ident(sequences[r.proteotype_id], sequences[c])
                if best is None or pid > best[0]:
                    best = (pid, diff, c)
            rows.append({"proteotype_id": r.proteotype_id, "family": fam,
                         "length": r.length, "n_cultivars": r.n_cultivars,
                         "nearest_reference_proteotype": best[2],
                         "identity_pct": round(best[0], 2),
                         "residue_differences": int(best[1])})
    return pd.DataFrame(rows)


def read_fasta(path) -> dict:
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0].split("__")[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


# --------------------------------------------------------------------------------------
# 4. prevalence at three granularities + bootstrap CI
# --------------------------------------------------------------------------------------
def _partition(prev: np.ndarray, n_cult: int) -> dict:
    soft = math.ceil(0.9 * n_cult)
    return {"core": int((prev == n_cult).sum()),
            "soft_core": int(((prev >= soft) & (prev < n_cult)).sum()),
            "shell": int(((prev >= 2) & (prev < soft)).sum()),
            "cloud": int((prev == 1).sum())}


def multi_granularity_partition(proteotypes: pd.DataFrame, members_g: pd.DataFrame,
                                cultivars: list) -> pd.DataFrame:
    """Core/soft-core/shell/cloud at transcript, gene and locus granularity.

    The same data give wildly different accessory fractions depending on how
    sequences are grouped, which is why a protein-sequence-level fraction cannot
    be compared with a DNA-level pan-genome core fraction.
    """
    n = len(cultivars)
    rows = []
    for fam in FAMILIES:
        # (a) transcript level, exact-identity proteotypes
        s = proteotypes[proteotypes["family"] == fam]
        p = _partition(s["n_cultivars"].to_numpy(), n)
        rows.append({"family": fam, "granularity": "transcript_exact_identity",
                     "clusters": len(s), **p})
        # (b) gene level: longest isoform per gene, clustered by its proteotype
        sub = members_g[members_g["family"] == fam]
        longest = sub.sort_values("length", ascending=False).groupby(["cultivar", "gene"]).head(1)
        tab = (longest.groupby(["proteotype_id", "cultivar"]).size().unstack(fill_value=0) > 0)
        tab = tab.reindex(columns=cultivars, fill_value=False)
        p = _partition(tab.sum(axis=1).to_numpy(), n)
        rows.append({"family": fam, "granularity": "gene_longest_isoform",
                     "clusters": int(tab.shape[0]), **p})
        # (c) locus level: loci matched by chromosome across cultivars
        per_cult = sub.groupby("cultivar")["gene"].nunique()
        chrom_sig = (sub.groupby(["cultivar", "chromosome"])["gene"].nunique()
                     .unstack(fill_value=0))
        invariant = bool((chrom_sig.nunique(axis=0) == 1).all())
        rows.append({"family": fam, "granularity": "gene_locus_chromosome_matched",
                     "clusters": int(per_cult.median()),
                     "core": int(per_cult.median()) if invariant else np.nan,
                     "soft_core": 0 if invariant else np.nan,
                     "shell": 0 if invariant else np.nan,
                     "cloud": 0 if invariant else np.nan})
    out = pd.DataFrame(rows)
    acc = out["shell"] + out["cloud"]
    out["accessory_pct"] = (100 * acc / out["clusters"]).round(1)
    # A locus row whose per-cultivar count is not invariant has no defined partition;
    # printing 0.0% there reads as "no accessory fraction", which is the opposite of
    # what the blank means.
    out["accessory_pct"] = out["accessory_pct"].astype(object).where(
        out["core"].notna(), "n/a")
    for c in ("core", "soft_core", "shell", "cloud"):
        out[c] = out[c].astype(object).where(out[c].notna(), "n/a")
    return out


def core_fraction_ci(proteotypes: pd.DataFrame, n_cult: int,
                     n_boot=10000, seed=1) -> pd.DataFrame:
    """Percentile bootstrap CI on the core fraction, resampling proteotypes."""
    rng = np.random.default_rng(seed)
    rows = []
    for fam in FAMILIES:
        prev = proteotypes.loc[proteotypes["family"] == fam, "n_cultivars"].to_numpy()
        if prev.size == 0:
            continue
        est = np.array([(prev[rng.integers(0, prev.size, prev.size)] == n_cult).mean()
                        for _ in range(n_boot)])
        lo, hi = np.percentile(est, [2.5, 97.5])
        rows.append({"family": fam, "n_proteotypes": int(prev.size),
                     "core_fraction_pct": round(100 * float((prev == n_cult).mean()), 1),
                     "ci95_low_pct": round(100 * float(lo), 1),
                     "ci95_high_pct": round(100 * float(hi), 1)})
    return pd.DataFrame(rows)




# --------------------------------------------------------------------------------------
# automatic detection of a divergent, chromosomally distinct clade
# --------------------------------------------------------------------------------------
def detect_divergent_clade(members_g: pd.DataFrame, family: str = "RCA",
                           nwk_path: str | None = None,
                           min_subgenomes: int = 3) -> dict:
    """Find proteotypes of `family` that sit on a minority homoeologous group.

    In wheat the canonical RCA genes are the tandem TaRca1/TaRca2 pair on
    homoeologous group 4, so group 4 holds twice as many loci as any other group.
    A second, complete homoeologous set on a different group is therefore a
    distinct lineage rather than an assembly artefact, and it is exactly the thing
    a candidate table should feature. Monophyly and support are then checked
    against the tree; they are reported, not assumed.
    """
    sub = members_g[(members_g["family"] == family)
                    & (members_g["homoeologous_group"] != "un")]
    if sub.empty:
        return {}
    loci = sub.groupby("homoeologous_group")["gene"].nunique().sort_values(ascending=False)
    if len(loci) < 2:
        return {}
    canonical = loci.index[0]
    out = {}
    for grp in loci.index[1:]:
        g = sub[sub["homoeologous_group"] == grp]
        subgen = set(g["subgenome"].unique()) - {"un"}
        if len(subgen) < min_subgenomes:
            continue
        pts = sorted(g["proteotype_id"].unique())
        if len(pts) < 2:
            continue
        out = {"family": family, "canonical_group": str(canonical),
               "clade_group": str(grp),
               "canonical_loci": int(loci.iloc[0]), "clade_loci": int(loci[grp]),
               "subgenomes": sorted(subgen), "members": pts, "n_members": len(pts)}
        break
    if not out:
        return {}
    if nwk_path and os.path.exists(nwk_path):
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            import tree_tools
            root = tree_tools.parse_newick(open(nwk_path).read())
            names = {l["name"] for l in tree_tools.leaves(root)}
            present = [p for p in out["members"] if p in names]
            if len(present) >= 2:
                mono, sup, extra, _ = tree_tools.clade_support(root, present)
                out.update(monophyletic=bool(mono),
                           support=(float(sup) if sup is not None else None),
                           extra_tips=sorted(extra), tested_members=present)
        except Exception:  # noqa: BLE001
            pass
    return out


# --------------------------------------------------------------------------------------
# 6. prioritised candidate table
# --------------------------------------------------------------------------------------
def candidate_table(proteotypes: pd.DataFrame, naming: pd.DataFrame,
                    members_g: pd.DataFrame, divergence: pd.DataFrame,
                    n_cult: int, extra_clade: list | None = None,
                    min_prevalence_fraction=0.8,
                    clade_label: str | None = None) -> pd.DataFrame:
    """Main-text candidate list: widely shared reference-absent proteotypes plus
    any explicitly named clade (e.g. the group-2 RCA-like clade)."""
    pt = proteotypes.copy()
    pt["novel"] = pt["novel_vs_reference"].astype(str).str.lower().isin(["true", "1"])
    thresh = math.ceil(min_prevalence_fraction * (n_cult - 1))
    widely = pt[(pt["novel"]) & (pt["n_cultivars"] >= thresh)]["proteotype_id"].tolist()
    clade = list(extra_clade or [])
    keep = sorted(set(widely) | set(clade))

    chrom = (members_g.groupby("proteotype_id")["chromosome"]
             .agg(lambda s: ";".join(sorted(set(s)))))
    t = (pt[pt["proteotype_id"].isin(keep)]
         .merge(naming[["proteotype_id", "uniprot_id", "representative_cultivar",
                        "composite_id"]], on="proteotype_id", how="left")
         .merge(divergence[["proteotype_id", "nearest_reference_proteotype",
                            "identity_pct", "residue_differences"]],
                on="proteotype_id", how="left"))
    t["chromosome"] = t["proteotype_id"].map(chrom)
    t["candidate_class"] = np.where(t["proteotype_id"].isin(clade),
                                    clade_label or "divergent clade",
                                    "reference-absent, widely shared")
    cols = ["proteotype_id", "composite_id", "candidate_class", "family", "chromosome",
            "length", "n_cultivars", "novel", "uniprot_id", "representative_cultivar",
            "nearest_reference_proteotype", "identity_pct", "residue_differences"]
    return t[cols].sort_values(["candidate_class", "family", "n_cultivars"],
                              ascending=[True, True, False])


# --------------------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------------------
def run_audit(inputs_dir, outdir, reference_cultivar="chinese_spring", group2_clade=None,
              query_table_hint=None):
    os.makedirs(outdir, exist_ok=True)

    def _find(fn):
        # reference_query_sets.csv is written to <project>/inputs/references/, which is
        # OUTSIDE the timestamped build directory, so searching the build alone misses it
        # and the Cpn60 classifier silently falls back to identity bands.
        up = os.path.abspath(os.path.join(inputs_dir, "..", ".."))
        for c in (os.path.join(inputs_dir, fn),
                  os.path.join(inputs_dir, "supplementary", fn),
                  os.path.join(inputs_dir, "references", fn),
                  os.path.join(up, "inputs", fn),
                  os.path.join(up, "inputs", "references", fn),
                  (os.path.join(query_table_hint, fn) if query_table_hint else "")):
            if c and os.path.exists(c):
                return c
        return None

    members = pd.read_csv(_find("members.csv"))
    proteotypes = pd.read_csv(_find("proteotypes.csv"))
    hits = pd.read_csv(_find("diamond_hits_all.csv"))
    naming_p = _find("proteotype_naming.csv")
    naming = pd.read_csv(naming_p) if naming_p else proteotypes.assign(
        uniprot_id="", representative_cultivar="", composite_id=proteotypes.proteotype_id)
    fasta = _find("PanProteome_catalog.annotated.faa") or _find("PanProteome_catalog.faa")

    members_g = resolve_gene_loci(members)
    cultivars = sorted(members_g["cultivar"].unique())
    n_cult = len(cultivars)
    out = {}

    members_g.to_csv(os.path.join(outdir, "gene_locus_table.csv"), index=False)
    out["locus_counts"] = locus_counts(members_g)
    out["locus_counts"].to_csv(os.path.join(outdir, "locus_counts_by_cultivar.csv"), index=False)
    out["cs_loci"] = cs_locus_table(members_g, reference_cultivar)
    out["cs_loci"].to_csv(os.path.join(outdir, "cs_locus_table.csv"), index=False)

    best, summary = specificity_audit(hits, members_g)
    best.to_csv(os.path.join(outdir, "specificity_per_model.csv"), index=False)
    summary.to_csv(os.path.join(outdir, "specificity_audit.csv"), index=False)
    out["specificity"] = summary

    if (members_g["family"] == "CPN60").any():
        cls = assign_cpn60_classes(
            best, members_g, hits=hits,
            query_table=(_find("reference_query_sets.csv") or ""))
        cls.to_csv(os.path.join(outdir, "cpn60_classes.csv"), index=False)
        out["cpn60_classes"] = cls

    if fasta:
        div = novel_divergence(proteotypes, read_fasta(fasta))
        div.to_csv(os.path.join(outdir, "novel_divergence.csv"), index=False)
        out["divergence"] = div
    else:
        div = pd.DataFrame(columns=["proteotype_id", "nearest_reference_proteotype",
                                    "identity_pct", "residue_differences"])
        sys.stderr.write("[audit] no proteotype FASTA found; divergence step skipped\n")

    out["granularity"] = multi_granularity_partition(proteotypes, members_g, cultivars)
    out["granularity"].to_csv(os.path.join(outdir, "granularity_partition.csv"), index=False)
    out["core_ci"] = core_fraction_ci(proteotypes, n_cult)
    out["core_ci"].to_csv(os.path.join(outdir, "core_fraction_ci.csv"), index=False)

    # If no clade was supplied, derive one. Proteotype ids change between runs, so a
    # hand-maintained list is guaranteed to go stale; detection from gene loci does not.
    clade_info = {}
    if not group2_clade:
        for fam in FAMILIES:
            nwk = _find(f"{fam}.nwk")
            ci = detect_divergent_clade(members_g, fam, nwk)
            if ci:
                clade_info = ci
                group2_clade = ci["members"]
                msg = (f"[audit] {fam}: detected a divergent clade on homoeologous group "
                       f"{ci['clade_group']} ({ci['clade_loci']} loci across subgenomes "
                       f"{'/'.join(ci['subgenomes'])}) against a canonical group "
                       f"{ci['canonical_group']} ({ci['canonical_loci']} loci); "
                       f"{ci['n_members']} proteotypes")
                if ci.get("support") is not None:
                    msg += (f"; monophyletic={ci.get('monophyletic')} "
                            f"support={ci['support']:.2f}")
                print(msg)
                break
    if clade_info:
        _flat = {k: (";".join(map(str, v)) if isinstance(v, (list, tuple)) else v)
                 for k, v in clade_info.items()}
        pd.DataFrame([_flat]).to_csv(
            os.path.join(outdir, "divergent_clade.csv"), index=False)
        out["clade_info"] = clade_info

    out["candidates"] = candidate_table(
        proteotypes, naming, members_g, div, n_cult, extra_clade=group2_clade,
        clade_label=(f"group-{clade_info['clade_group']} "
                     f"{clade_info['family']}-like clade" if clade_info else None))
    out["candidates"].to_csv(os.path.join(outdir, "candidate_table.csv"), index=False)

    print(f"[audit] {n_cult} cultivars")
    for fam in FAMILIES:
        s = members_g[members_g["family"] == fam]
        if s.empty:
            continue
        per = s.groupby("cultivar")["gene"].nunique()
        print(f"[audit] {fam}: {s.source_id.nunique()} protein models -> "
              f"{s.gene.nunique()} gene loci "
              f"({per.min()}-{per.max()} per cultivar"
              f"{', INVARIANT' if per.nunique() == 1 else ''})")
    print(f"[audit] outputs -> {os.path.abspath(outdir)}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", required=True, help="pan_proteome_build_* directory")
    ap.add_argument("--outdir", default=None, help="default <inputs>/audit")
    ap.add_argument("--reference", default="chinese_spring")
    ap.add_argument("--query-table-dir", default=None,
                    help="directory holding reference_query_sets.csv, if not auto-found")
    ap.add_argument("--group2-clade", default="",
                    help="comma-separated proteotype ids to force into the candidate table")
    a = ap.parse_args(argv)
    clade = [x.strip() for x in a.group2_clade.split(",") if x.strip()]
    run_audit(a.inputs, a.outdir or os.path.join(a.inputs, "audit"),
              reference_cultivar=a.reference, group2_clade=clade or None,
              query_table_hint=a.query_table_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
