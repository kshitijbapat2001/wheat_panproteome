#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

import pandas as pd

DEFAULT_ANCHORS = {
    "TaRca1": {"A": "TraesCS4A02G177600", "B": "TraesCS4B02G140200", "D": "TraesCS4D02G134900"},
    "TaRca2": {"A": "TraesCS4A02G177500", "B": "TraesCS4B02G140300", "D": "TraesCS4D02G135000"},
}
DEFAULT_HRE = {("TaRca1", "A"): "reported", ("TaRca1", "B"): "reported",
               ("TaRca1", "D"): "reported", ("TaRca2", "A"): "reported",
               ("TaRca2", "B"): "not detected", ("TaRca2", "D"): "not detected"}

ISO_LABEL = {"RCA1b": "RCA1β", "RCA2a": "RCA2α", "RCA2b": "RCA2β",
             "RCAL": "RCA-like"}
ISO_ORDER = ["RCA1b", "RCA2a", "RCA2b", "RCAL"]
SUBS = ["A", "B", "D"]


def _find(root: str, name: str):
    for c in (os.path.join(root, name), os.path.join(root, "audit", name),
              os.path.join(root, "domains", name),
              os.path.join(root, "final", "supplementary", "data", name)):
        if os.path.exists(c):
            return c
    hits = glob.glob(os.path.join(root, "**", name), recursive=True)
    return sorted(hits, key=os.path.getmtime, reverse=True)[0] if hits else None


def read_fasta(path: str) -> dict:
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0].split("__")[0]
                buf = []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def _truthy(x) -> bool:
    return str(x).strip().lower() in ("true", "1", "yes")


def make_aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    for new, old in (("end_insertion_score", "target_end_gap_score"),
                     ("end_deletion_score", "query_end_gap_score")):
        try:
            setattr(al, new, 0.0)
        except AttributeError:
            setattr(al, old, 0.0)
    return al


def make_local_aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.mode = "local"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    return al


def local_map(al, seg: str, qry: str) -> dict:
    alns = al.align(seg, qry)
    try:
        aln = alns[0]
    except IndexError:
        return {"pairs": {}, "identity": 0.0}
    pairs = {}
    for (t0, t1), (q0, q1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(int(t1) - int(t0)):
            pairs[int(t0) + k] = int(q0) + k
    same = sum(1 for t, q in pairs.items() if seg[t] == qry[q])
    return {"pairs": pairs, "identity": 100.0 * same / len(pairs) if pairs else 0.0}


def pair_map(al, ref: str, qry: str) -> dict:
    aln = al.align(ref, qry)[0]
    tb, qb = aln.aligned
    pairs = {}
    tb = [(int(a), int(b)) for a, b in tb]
    qb = [(int(a), int(b)) for a, b in qb]
    for (t0, t1), (q0, q1) in zip(tb, qb):
        for k in range(t1 - t0):
            pairs[t0 + k] = q0 + k
    if not pairs:
        return {"pairs": {}, "identity": 0.0, "n_aligned": 0, "coverage": 0.0,
                "n_missing": len(ref), "c_missing": 0, "del_blocks": [],
                "ins_blocks": [], "subs": []}
    ref_first, ref_last = min(pairs), max(pairs)
    same = sum(1 for t, q in pairs.items() if ref[t] == qry[q])
    subs = [(t, ref[t], qry[q]) for t, q in sorted(pairs.items()) if ref[t] != qry[q]]
    dels, run = [], None
    for t in range(ref_first, ref_last + 1):
        if t not in pairs:
            run = [t, t] if run is None else [run[0], t]
        elif run is not None:
            dels.append(tuple(run)); run = None
    if run is not None:
        dels.append(tuple(run))
    ins = []
    for (a, b) in zip(qb[:-1], qb[1:]):
        if b[0] > a[1]:
            ins.append((a[1], b[0] - a[1]))
    return {"pairs": pairs, "identity": 100.0 * same / len(pairs),
            "n_aligned": len(pairs), "coverage": len(pairs) / len(ref),
            "n_missing": ref_first, "c_missing": len(ref) - 1 - ref_last,
            "del_blocks": dels, "ins_blocks": ins, "subs": subs,
            "q_first": min(pairs.values()), "q_last": max(pairs.values())}


def assign(pt: pd.DataFrame, mem: pd.DataFrame, seqs: dict, naming: pd.DataFrame | None,
           anchors: dict, hre: dict, cte_min: int = 10, full_cov: float = 0.97,
           max_len_dev: float = 0.03, max_indel: int = 4):
    al = make_aligner()
    rca = pt[pt["family"] == "RCA"].copy()
    rca["novel"] = rca["novel_vs_reference"].map(_truthy)
    rca["n_cultivars"] = rca["n_cultivars"].astype(int)
    m = mem[mem["family"] == "RCA"].copy()
    m["homoeologous_group"] = m["homoeologous_group"].astype(str)

    anchor_of = {g: (gene, sub) for gene, d in anchors.items() for sub, g in d.items()}
    loci = m.groupby("proteotype_id").agg(
        genes=("gene", lambda s: sorted(set(s))),
        chromosomes=("chromosome", lambda s: sorted(set(map(str, s)))),
        subgenomes=("subgenome", lambda s: sorted(set(map(str, s)))),
        groups=("homoeologous_group", lambda s: sorted(set(s)))).to_dict("index")

    info = {}
    for r in rca.itertuples():
        L = loci.get(r.proteotype_id, {})
        info[r.proteotype_id] = {"seq": seqs.get(r.proteotype_id, ""), "n": r.n_cultivars,
                                 "novel": r.novel, "cultivars": r.cultivars,
                                 "genes": L.get("genes", []),
                                 "chromosomes": L.get("chromosomes", []),
                                 "subgenomes": L.get("subgenomes", []),
                                 "groups": L.get("groups", [])}
    missing = [p for p, v in info.items() if not v["seq"]]
    if missing:
        raise SystemExit(f"[rca_isoforms] no sequence for {len(missing)} proteotypes, "
                         f"e.g. {missing[:3]}; check the catalogue FASTA")

    anchor_pts = {"TaRca1": [], "TaRca2": []}
    for p, v in info.items():
        for g in v["genes"]:
            if g in anchor_of:
                anchor_pts[anchor_of[g][0]].append(p)
    for k in anchor_pts:
        anchor_pts[k] = sorted(set(anchor_pts[k]))
    if not anchor_pts["TaRca1"] or not anchor_pts["TaRca2"]:
        raise SystemExit("[rca_isoforms] Chinese Spring anchor loci not found in members.csv; "
                         "check rca_isoforms.anchors in config.yaml")

    def best_identity(p, pool):
        best = (0.0, None)
        for c in pool:
            if c == p:
                continue
            a = pair_map(al, info[c]["seq"], info[p]["seq"])
            short = min(len(info[c]["seq"]), len(info[p]["seq"]))
            if a["n_aligned"] >= 0.5 * short and a["identity"] > best[0]:
                best = (a["identity"], c)
        return best

    for p, v in info.items():
        direct = sorted({anchor_of[g][0] for g in v["genes"] if g in anchor_of})
        v["id_rca1"], v["id_rca2"] = best_identity(p, anchor_pts["TaRca1"])[0], \
                                     best_identity(p, anchor_pts["TaRca2"])[0]
        if v["groups"] == ["2"]:
            v["gene"], v["gene_evidence"] = "RCA-like (group 2)", "homoeologous group 2 locus"
        elif "2" in v["groups"]:
            v["gene"], v["gene_evidence"] = "unresolved", "loci on groups 2 and 4"
        elif len(direct) == 1:
            v["gene"] = direct[0]
            v["gene_evidence"] = "Chinese Spring anchor locus (Caruana et al. 2022)"
        elif len(direct) > 1:
            v["gene"], v["gene_evidence"] = "unresolved", "maps to both anchor genes"
        else:
            v["gene"] = "TaRca1" if v["id_rca1"] > v["id_rca2"] else "TaRca2"
            v["gene_evidence"] = (f"sequence: {max(v['id_rca1'], v['id_rca2']):.1f}% vs "
                                  f"{min(v['id_rca1'], v['id_rca2']):.1f}% identity to the "
                                  "Chinese Spring TaRca1/TaRca2 proteotypes")

    gene_by_locus = {}
    for p, v in info.items():
        for g in v["genes"]:
            gene_by_locus.setdefault(g, set()).add(v["gene"])
    inconsistent = {g: s for g, s in gene_by_locus.items() if len(s) > 1}

    def last_cys(s, w=40):
        return s[-w:].count("C")

    def prevalence_key(p):
        v = info[p]
        return (-v["n"], v["novel"], -len(v["seq"]), p)

    rca2 = [p for p, v in info.items() if v["gene"] == "TaRca2"]
    alpha_ref = sorted([p for p in rca2 if last_cys(info[p]["seq"]) >= 2], key=prevalence_key)[0]
    a_seq = info[alpha_ref]["seq"]
    beta_ref = sorted([p for p in rca2 if last_cys(info[p]["seq"]) == 0
                       and len(info[p]["seq"]) >= 0.85 * len(a_seq)], key=prevalence_key)[0]
    b_seq = info[beta_ref]["seq"]
    ab = pair_map(al, a_seq, b_seq)
    cte_start = max(ab["pairs"]) + 1
    cte_cys = [cte_start + i for i, c in enumerate(a_seq[cte_start:]) if c == "C"]

    loc = make_local_aligner()
    cte_seg = a_seq[cte_start:]
    for p, v in info.items():
        lm = local_map(loc, cte_seg, v["seq"])
        has_cte = len(lm["pairs"]) >= cte_min and lm["identity"] >= 60.0
        v["cte_aligned"] = len(lm["pairs"]) if has_cte else 0
        v["cte_q"] = sorted(lm["pairs"].values()) if has_cte else []
        v["cte_t0"] = min(lm["pairs"]) if has_cte else None
        a = {"pairs": {cte_start + t: q for t, q in lm["pairs"].items()}} if has_cte else {"pairs": {}}
        if v["gene"] == "TaRca2":
            v["iso"] = "RCA2a" if has_cte else "RCA2b"
        elif v["gene"] == "TaRca1":
            v["iso"] = "RCA1b"
            if has_cte:
                v["gene_evidence"] += "; WARNING: CTE-like sequence in a TaRca1 model"
        elif v["gene"].startswith("RCA-like"):
            v["iso"] = "RCAL"
        else:
            v["iso"] = "unresolved"
        v["k428"] = v["seq"][a["pairs"][cte_start]] if cte_start in a["pairs"] else ""
        v["cte_cys"] = "/".join(v["seq"][a["pairs"][c]] if c in a["pairs"] else "-"
                                for c in cte_cys) if has_cte else ""

    mm = re.search(r"EAAD([MIL])IKK", b_seq)
    if not mm:
        raise SystemExit("[rca_isoforms] could not locate the Degen et al. (2020) switch "
                         "site motif in the canonical RCA2beta")
    m159 = mm.start(1)
    for p, v in info.items():
        b = pair_map(al, b_seq, v["seq"])
        v["m159"] = v["seq"][b["pairs"][m159]] if m159 in b["pairs"] else "-"
        v["id_rca2b"] = b["identity"]

    canon = {}
    for iso in ISO_ORDER:
        members = [p for p, v in info.items() if v["iso"] == iso]
        if members:
            canon[iso] = sorted(members, key=prevalence_key)[0]

    def describe(p, ref_id):
        ref, q = info[ref_id]["seq"], info[p]["seq"]
        a = pair_map(al, ref, q)
        big_del = [(s, e) for s, e in a["del_blocks"] if e - s + 1 > max_indel]
        big_ins = [(qp, n) for qp, n in a["ins_blocks"] if n > max_indel]
        n_extra = a.get("q_first", 0)
        c_extra = len(q) - 1 - a.get("q_last", len(q) - 1)
        notes = []
        if a["n_missing"] > max_indel:
            notes.append(f"N-terminal {a['n_missing']} residues absent")
        if n_extra > max_indel:
            notes.append(f"alternative N-terminus (+{n_extra})")
        if big_del:
            notes.append("internal deletion " + ", ".join(
                f"{s + 1}-{e + 1}" for s, e in big_del))
        if big_ins:
            notes.append("internal insertion " + ", ".join(f"{n} aa" for _, n in big_ins))
        if a["c_missing"] > max_indel:
            notes.append(f"C-terminal {a['c_missing']} residues absent")
        if c_extra > max_indel:
            notes.append(f"alternative C-terminus (+{c_extra})")
        full = (not notes and a["coverage"] >= full_cov
                and abs(len(q) - len(ref)) <= max_len_dev * len(ref))
        return a, full, "; ".join(notes)

    for p, v in info.items():
        if v["iso"] not in canon:
            v.update(full=False, model="", cov=float("nan"))
            continue
        a, full, notes = describe(p, canon[v["iso"]])
        if v["iso"] == "RCA2a" and v["cte_q"] and "C-terminal" in notes:
            ref = info[canon[v["iso"]]]["seq"]
            q0, t0 = v["cte_q"][0], cte_start + v["cte_t0"]
            body = pair_map(al, ref[:t0], v["seq"][:q0])
            last = max(body["pairs"]) if body["pairs"] else -1
            parts = [n for n in notes.split("; ") if not n.startswith("C-terminal")]
            if t0 - last - 1 > max_indel:
                parts.append(f"internal deletion {last + 2}-{t0} (CTE tail retained)")
            notes = "; ".join(parts)
        v["full"], v["model"], v["cov"] = full, ("full-length" if full else
                                                 "variant model: " + (notes or "partial")), a["coverage"]

    def sub_label(v):
        subs = sorted({s for s in v["subgenomes"] if s in SUBS})
        return subs

    for p, v in info.items():
        subs = sub_label(v)
        grp = "4" if v["iso"] in ("RCA1b", "RCA2a", "RCA2b") else "2"
        chrs = "/".join(f"{grp}{s}" for s in subs) if subs else "un"
        v["label"] = f"{ISO_LABEL.get(v['iso'], v['iso'])}-{chrs}" if v["iso"] != "unresolved" \
            else "unresolved"
        v["cells"] = [(v["iso"], s) for s in subs]

    cell_canon = {}
    for iso in ISO_ORDER:
        for s in SUBS:
            full = [p for p, v in info.items() if (iso, s) in v["cells"] and v["full"]]
            if full:
                cell_canon[(iso, s)] = sorted(full, key=prevalence_key)[0]

    for p, v in info.items():
        refs = [cell_canon[c] for c in v["cells"] if c in cell_canon]
        if not refs:
            v.update(ref_cell="", n_sub="", subs="", id_cell="")
            continue
        ref_id = refs[0]
        a = pair_map(al, info[ref_id]["seq"], v["seq"])
        v["ref_cell"] = ref_id
        v["n_sub"] = len(a["subs"])
        v["subs"] = "; ".join(f"{r}{t + 1}{q}" for t, r, q in a["subs"])
        v["id_cell"] = round(a["identity"], 2)

    comp = {}
    if naming is not None and "composite_id" in naming.columns:
        comp = dict(zip(naming["proteotype_id"], naming["composite_id"]))
    rows = []
    for p in sorted(info, key=lambda x: int(x.split("PT")[-1])):
        v = info[p]
        rows.append({
            "proteotype_id": p,
            "composite_id": comp.get(p, ""),
            "established_isoform": v["label"],
            "gene": v["gene"],
            "isoform": ISO_LABEL.get(v["iso"], v["iso"]),
            "chromosomes": ";".join(v["chromosomes"]),
            "gene_loci": ";".join(v["genes"]),
            "assignment_evidence": v["gene_evidence"],
            "identity_to_TaRca1_pct": round(v["id_rca1"], 2),
            "identity_to_TaRca2_pct": round(v["id_rca2"], 2),
            "cte_residues_aligned": int(v["cte_aligned"]),
            "length": len(v["seq"]),
            "model_status": v["model"],
            "coverage_of_canonical_pct": round(100 * v["cov"], 1) if v["cov"] == v["cov"] else "",
            "canonical_proteotype": v["ref_cell"],
            "identity_to_canonical_pct": v["id_cell"],
            "substitutions_vs_canonical": v["n_sub"],
            "substitution_list": v["subs"],
            "switch_site_M159I": v["m159"],
            "site_K428": v["k428"] if v["iso"] == "RCA2a" else "",
            "cte_cysteine_pair": v["cte_cys"] if v["iso"] == "RCA2a" else "",
            "identity_to_RCA2beta_pct": round(v["id_rca2b"], 2),
            "n_cultivars": v["n"],
            "cultivars": v["cultivars"],
            "novel_vs_reference": v["novel"],
        })
    per = pd.DataFrame(rows)[[
        "proteotype_id", "composite_id", "established_isoform", "gene", "isoform",
        "chromosomes", "gene_loci", "length", "n_cultivars", "cultivars",
        "novel_vs_reference", "model_status", "coverage_of_canonical_pct",
        "canonical_proteotype", "identity_to_canonical_pct", "substitutions_vs_canonical",
        "substitution_list", "switch_site_M159I", "site_K428", "cte_cysteine_pair",
        "cte_residues_aligned", "identity_to_TaRca1_pct", "identity_to_TaRca2_pct",
        "identity_to_RCA2beta_pct", "assignment_evidence"]]

    cs_locus = {}
    for gene, d in anchors.items():
        for s, g in d.items():
            cs_locus[(("RCA1b",) if gene == "TaRca1" else ("RCA2a", "RCA2b")), s] = g
    grp2_cs = {}
    for p, v in info.items():
        if v["iso"] == "RCAL":
            for g, c in zip(v["genes"], v["chromosomes"] * len(v["genes"])):
                if g.startswith("TraesCS"):
                    grp2_cs.setdefault(g[7:9], g)

    summ = []
    for iso in ISO_ORDER:
        for s in SUBS:
            ps = [p for p, v in info.items() if (iso, s) in v["cells"]]
            if iso == "RCAL" and not ps:
                continue
            full = [p for p in ps if info[p]["full"]]
            fcult = set()
            for p in full:
                fcult |= set(str(info[p]["cultivars"]).split(";"))
            if iso == "RCA1b":
                locus, gene = anchors["TaRca1"][s], "TaRca1"
            elif iso in ("RCA2a", "RCA2b"):
                locus, gene = anchors["TaRca2"][s], "TaRca2"
            else:
                locus, gene = grp2_cs.get(f"2{s}", ""), "RCA-like"
            princ = cell_canon.get((iso, s), "")
            summ.append({
                "established_isoform": f"{ISO_LABEL[iso]}-{'2' if iso == 'RCAL' else '4'}{s}",
                "gene": gene,
                "chinese_spring_locus": locus,
                "proteoforms": len(ps),
                "full_length": len(full),
                "reference_absent": sum(1 for p in ps if info[p]["novel"]),
                "principal_full_length_proteoform": princ,
                "principal_n_cultivars": info[princ]["n"] if princ else 0,
                "cultivars_with_full_length_form": len(fcult - {""}),
                "chinese_spring_full_length": "chinese_spring" in fcult,
                "promoter_heat_responsive_element": hre.get((gene, s), "not examined"),
                "proteoform_ids": ";".join(sorted(ps, key=lambda x: int(x.split("PT")[-1]))),
            })
    summ = pd.DataFrame(summ)
    ref = {"canonical": {ISO_LABEL[k]: v for k, v in canon.items()},
           "canonical_by_cell": {f"{ISO_LABEL[k[0]]}-{k[1]}": v for k, v in cell_canon.items()},
           "alpha_reference": alpha_ref, "beta_reference": beta_ref,
           "cte_start_1based": cte_start + 1, "cte_length": len(a_seq) - cte_start,
           "cte_cysteines_1based": [c + 1 for c in cte_cys],
           "switch_site_1based_on_beta_reference": m159 + 1,
           "anchors": anchors, "inconsistent_loci": {k: sorted(v) for k, v in inconsistent.items()}}
    return per, summ, ref


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", required=True, help="pan_proteome_build_* directory")
    ap.add_argument("--outdir", default=None, help="default <inputs>/audit")
    ap.add_argument("--config", default="config.yaml")
    a = ap.parse_args(argv)
    outdir = a.outdir or os.path.join(a.inputs, "audit")
    os.makedirs(outdir, exist_ok=True)

    cfg = {}
    if os.path.exists(a.config):
        import yaml
        cfg = (yaml.safe_load(open(a.config)) or {}).get("rca_isoforms", {}) or {}
    anchors = cfg.get("anchors") or DEFAULT_ANCHORS
    hre = dict(DEFAULT_HRE)
    for k, v in (cfg.get("heat_responsive_elements") or {}).items():
        gene, sub = k.split("-")
        hre[(gene, sub)] = v

    pt_p, mem_p = _find(a.inputs, "proteotypes.csv"), _find(a.inputs, "members.csv")
    fa_p = (_find(a.inputs, "PanProteome_catalog.faa")
            or _find(a.inputs, "PanProteome_catalog.annotated.faa"))
    if not (pt_p and mem_p and fa_p):
        raise SystemExit(f"[rca_isoforms] missing inputs: proteotypes={pt_p} members={mem_p} "
                         f"catalogue={fa_p}")
    nm_p = _find(a.inputs, "proteotype_naming.csv")
    per, summ, ref = assign(pd.read_csv(pt_p), pd.read_csv(mem_p), read_fasta(fa_p),
                            pd.read_csv(nm_p) if nm_p else None, anchors, hre,
                            cte_min=int(cfg.get("cte_min_residues", 10)),
                            full_cov=float(cfg.get("full_length_min_coverage", 0.97)))
    per.to_csv(os.path.join(outdir, "rca_isoform_assignment.csv"), index=False)
    summ.to_csv(os.path.join(outdir, "rca_isoform_summary.csv"), index=False)
    with open(os.path.join(outdir, "rca_isoform_reference.json"), "w") as fh:
        json.dump(ref, fh, indent=2)
    n_full = int(per["model_status"].eq("full-length").sum())
    print(f"[rca_isoforms] {len(per)} RCA proteotypes -> "
          f"{per['established_isoform'].nunique()} isoform labels; {n_full} full-length")
    print(f"[rca_isoforms] canonical RCA2alpha={ref['alpha_reference']} "
          f"RCA2beta={ref['beta_reference']}; CTE starts at residue {ref['cte_start_1based']}")
    if ref["inconsistent_loci"]:
        print(f"[rca_isoforms] WARNING loci with mixed gene calls: {ref['inconsistent_loci']}")
    print(f"[rca_isoforms] outputs -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
