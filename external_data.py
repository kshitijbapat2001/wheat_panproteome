#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""external_data.py - enrich the pan-proteome with primary database records"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

UNIPROT = "https://rest.uniprot.org/uniprotkb"
INTERPRO = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = "wheat-panproteome/3.0 (+https://github.com/)"

_LAST_CALL = [0.0]
_MIN_INTERVAL = 0.34


def _throttle():
    dt = time.time() - _LAST_CALL[0]
    if dt < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - dt)
    _LAST_CALL[0] = time.time()


def _get(url, timeout=60, retries=3, log=None):
    last = None
    for attempt in range(1, retries + 1):
        try:
            _throttle()
            req = urllib.request.Request(url, headers={"User-Agent": UA,
                                                       "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            last = e
            if log:
                log.debug("GET failed (%d/%d) %s: %s", attempt, retries, url, e)
            time.sleep(1.5 * attempt)
    if log:
        log.warning("giving up on %s (%s)", url, last)
    return None


def _cache_path(cache_dir, kind, key):
    os.makedirs(os.path.join(cache_dir, kind), exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in str(key))
    return os.path.join(cache_dir, kind, safe)


def _cached(cache_dir, kind, key, fetch, log=None):
    p = _cache_path(cache_dir, kind, key)
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return open(p, encoding="utf-8").read()
    txt = fetch()
    if txt:
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(txt)
    return txt


def uniprot_record(acc, cache_dir, log=None) -> dict:
    """Full UniProtKB JSON record for one accession (cached)."""
    txt = _cached(cache_dir, "uniprot_json", acc + ".json",
                  lambda: _get(f"{UNIPROT}/{acc}.json", log=log), log)
    if not txt:
        return {}
    try:
        return json.loads(txt)
    except Exception:
        return {}


def uniprot_fasta(acc, cache_dir, log=None) -> str:
    return _cached(cache_dir, "uniprot_fasta", acc + ".fasta",
                   lambda: _get(f"{UNIPROT}/{acc}.fasta", log=log), log) or ""


def summarise_uniprot(rec: dict) -> dict:
    """Flatten a UniProtKB record into the fields a reviewer actually asks for."""
    if not rec:
        return {}
    d = rec.get("proteinDescription", {}) or {}
    name = ""
    for key in ("recommendedName", "submissionNames"):
        v = d.get(key)
        if isinstance(v, dict):
            name = v.get("fullName", {}).get("value", "")
        elif isinstance(v, list) and v:
            name = v[0].get("fullName", {}).get("value", "")
        if name:
            break
    genes = ";".join(g.get("geneName", {}).get("value", "")
                     for g in rec.get("genes", []) if g.get("geneName"))
    loc = []
    for c in rec.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for s in c.get("subcellularLocations", []):
                v = s.get("location", {}).get("value")
                if v:
                    loc.append(v)
    pfam, interpro, panther = [], [], []
    for x in rec.get("uniProtKBCrossReferences", []):
        db, xid = x.get("database"), x.get("id")
        props = {p.get("key"): p.get("value") for p in x.get("properties", [])}
        label = f"{xid}({props.get('EntryName', '')})" if props.get("EntryName") else xid
        if db == "Pfam":
            pfam.append(label)
        elif db == "InterPro":
            interpro.append(label)
        elif db == "PANTHER":
            panther.append(xid)
    return {"protein_name": name, "gene_names": genes,
            "subcellular_location": ";".join(dict.fromkeys(loc)),
            "pfam": ";".join(pfam), "interpro": ";".join(interpro),
            "panther": ";".join(dict.fromkeys(panther)),
            "sequence_length": rec.get("sequence", {}).get("length"),
            "reviewed": str(rec.get("entryType", "")).lower().startswith("uniprotkb reviewed")}


def interpro_entries(acc, cache_dir, log=None) -> str:
    """InterPro entries matched to a UniProt accession (cached, best effort)."""
    txt = _cached(cache_dir, "interpro", acc + ".json",
                  lambda: _get(f"{INTERPRO}/{acc}?page_size=100", log=log), log)
    if not txt:
        return ""
    try:
        js = json.loads(txt)
    except Exception:
        return ""
    out = []
    for r in js.get("results", []):
        md = r.get("metadata", {})
        if md.get("accession"):
            out.append(f"{md['accession']}({md.get('name', '')})")
    return ";".join(out)


def ncbi_esummary(db, term, cache_dir, email="", api_key="", log=None) -> dict:
    q = {"db": db, "term": term, "retmode": "json", "retmax": "1"}
    if email:
        q["email"] = email
    if api_key:
        q["api_key"] = api_key
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(q)
    txt = _cached(cache_dir, "ncbi_esearch", f"{db}_{term}.json", lambda: _get(url, log=log), log)
    if not txt:
        return {}
    try:
        ids = json.loads(txt).get("esearchresult", {}).get("idlist", [])
    except Exception:
        return {}
    if not ids:
        return {}
    q2 = {"db": db, "id": ids[0], "retmode": "json"}
    if email:
        q2["email"] = email
    if api_key:
        q2["api_key"] = api_key
    url2 = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(q2)
    txt2 = _cached(cache_dir, "ncbi_esummary", f"{db}_{ids[0]}.json",
                   lambda: _get(url2, log=log), log)
    try:
        return json.loads(txt2).get("result", {}).get(ids[0], {})
    except Exception:
        return {}


def annotate_naming_table(naming_csv, cache_dir, cfg, log=None):
    """Add UniProt/InterPro annotation columns to proteotype_naming.csv."""
    if not os.path.exists(naming_csv):
        if log:
            log.warning("[external] %s not found; annotation skipped", naming_csv)
        return None
    nm = pd.read_csv(naming_csv)
    accs = [a for a in nm.get("uniprot_id", pd.Series(dtype=str)).dropna().unique()
            if isinstance(a, str) and a]
    cap = int(cfg.get("external", {}).get("max_accessions", 400))
    if len(accs) > cap:
        if log:
            log.warning("[external] %d accessions exceeds max_accessions=%d; truncating",
                        len(accs), cap)
        accs = accs[:cap]
    if log:
        log.info("[external] annotating %d distinct UniProt accessions "
                 "(cache: %s)", len(accs), cache_dir)
    rows = []
    for i, acc in enumerate(accs, 1):
        rec = summarise_uniprot(uniprot_record(acc, cache_dir, log))
        rec["uniprot_id"] = acc
        rec["interpro_full"] = interpro_entries(acc, cache_dir, log)
        rows.append(rec)
        if log and i % 25 == 0:
            log.info("[external]   %d/%d", i, len(accs))
    ann = pd.DataFrame(rows)
    if ann.empty:
        return None
    out = nm.merge(ann, on="uniprot_id", how="left")
    return out


def append_outgroups(family_fasta, accessions, cache_dir, cfg, log=None):
    """Append outgroup sequences to a family FASTA before alignment."""
    if not accessions:
        return []
    have = set()
    if os.path.exists(family_fasta):
        with open(family_fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    have.add(line[1:].split()[0])
    added = []
    with open(family_fasta, "a") as fo:
        for acc in accessions:
            rid = f"OUTGROUP_{acc}"
            if rid in have:
                added.append(rid)
                continue
            fa = uniprot_fasta(acc, cache_dir, log)
            if not fa or not fa.startswith(">"):
                if log:
                    log.warning("[external] outgroup %s could not be fetched", acc)
                continue
            seq = "".join(l.strip() for l in fa.splitlines()[1:])
            if not seq:
                continue
            rec = summarise_uniprot(uniprot_record(acc, cache_dir, log))
            fo.write(f">{rid} outgroup uniprot={acc} "
                     f"name=\"{rec.get('protein_name', '')}\"\n")
            for i in range(0, len(seq), 60):
                fo.write(seq[i:i + 60] + "\n")
            added.append(rid)
    return added


def warm_cache(cfg, cache_dir, log=None):
    """Pre-fetch everything the compute nodes will need, on a machine with internet."""
    os.makedirs(cache_dir, exist_ok=True)
    ext = cfg.get("external", {})
    n = 0
    for fam, accs in (ext.get("outgroups") or {}).items():
        for acc in accs:
            uniprot_fasta(acc, cache_dir, log)
            uniprot_record(acc, cache_dir, log)
            n += 1
    print(f"[external] warmed {n} outgroup record(s) into {cache_dir}")
    print("[external] proteotype annotations are cached during the run; re-run "
          "--warm-cache after a build to make the cache complete for offline reuse.")
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--warm-cache", action="store_true")
    ap.add_argument("--annotate", default=None, help="path to proteotype_naming.csv")
    ap.add_argument("--out", default="proteotype_annotation.csv")
    a = ap.parse_args(argv)

    cfg = {}
    if os.path.exists(a.config):
        import yaml
        cfg = yaml.safe_load(open(a.config)) or {}
    cache = (a.cache_dir or cfg.get("external", {}).get("cache_dir")
             or os.path.join(cfg.get("project", {}).get("inputs_dir", "./inputs"),
                             "external_cache"))
    os.makedirs(cache, exist_ok=True)

    if a.warm_cache:
        warm_cache(cfg, cache)
    if a.annotate:
        ann = annotate_naming_table(a.annotate, cache, cfg)
        if ann is not None:
            ann.to_csv(a.out, index=False)
            print(f"[external] wrote {a.out} ({len(ann)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
