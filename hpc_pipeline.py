#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime
import gzip
import json
import logging
import math
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import textwrap
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 600, "savefig.dpi": 600})

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

def _panel_letter(ax, letter, dx=-52, dy=16):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy),
                textcoords="offset points", fontsize=14, fontweight="bold",
                va="bottom", ha="left")
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap

try:
    import seaborn as sns
    _HAVE_SNS = True
except Exception:
    _HAVE_SNS = False

try:
    from scipy import stats as scipy_stats
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

try:
    from Bio import SeqIO, AlignIO, Phylo
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        _HAVE_PROTPARAM = True
    except Exception:
        _HAVE_PROTPARAM = False
    _HAVE_BIO = True
except Exception:
    _HAVE_BIO = False
    _HAVE_PROTPARAM = False


STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

AA_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
WATER = 18.01524

PK = {"Nterm": 8.6, "Cterm": 3.6, "C": 8.5, "D": 3.9, "E": 4.1,
      "H": 6.5, "K": 10.8, "R": 12.5, "Y": 10.1}

HYDROPHOBIC = set("AVILMFWC")
FAMILIES = ("RCA", "CPN60")
FAM_DISPLAY = {"RCA": "RCA", "CPN60": "Cpn60"}
FIG_BLUE, FIG_RED, FIG_ABSENT, FIG_GREY = "#3B6FB6", "#C24B4B", "#F2F3F5", "#B8BEC7"
FIG_AMBER = "#E4A340"
CULT_LABELS = {
    "chinese_spring": "Chinese Spring", "arinalrfor": "ArinaLrFor", "jagger": "Jagger",
    "julius": "Julius", "lancer": "Lancer", "landmark": "Landmark", "mace": "Mace",
    "norin61": "Norin61", "stanley": "CDC Stanley", "mattis": "SY Mattis",
}
FAM_DISP = {"RCA": "RCA", "CPN60": "Cpn60"}


def get_timestamp() -> str:
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def sanitize_name(name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9]", "_", str(name))
    return re.sub(r"_+", "_", clean).strip("_") or "unnamed"


def setup_logging(out_dir: str) -> logging.Logger:
    os.makedirs(out_dir, exist_ok=True)
    log = logging.getLogger("panproteome")
    log.setLevel(logging.DEBUG)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)
    log.addHandler(sh)

    fh = logging.FileHandler(os.path.join(out_dir, "pipeline.log"))
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    log.addHandler(fh)
    return log


def require_tool(name: str, log: logging.Logger) -> str:
    path = shutil.which(name)
    if not path:
        raise RuntimeError(
            f"Required tool '{name}' not found on PATH. Activate the conda environment "
            f"(environment.yml) before running, e.g. `conda activate wheat_panproteome`."
        )
    log.debug("tool %-10s -> %s", name, path)
    return path


def run_cmd(cmd, log, cwd=None, env=None, check=True, quiet=False):
    printable = " ".join(str(c) for c in cmd)
    log.debug("RUN: %s", printable)
    proc = subprocess.run(
        [str(c) for c in cmd], cwd=cwd, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if proc.stdout and not quiet:
        for line in proc.stdout.splitlines():
            log.debug("    %s", line)
    if check and proc.returncode != 0:
        tail = "\n".join((proc.stdout or "").splitlines()[-25:])
        log.error("Command failed (rc=%s): %s\n%s", proc.returncode, printable, tail)
        raise RuntimeError(f"Command failed: {cmd[0]} (rc={proc.returncode})")
    return proc


def tool_version(name: str, args, log) -> str:
    try:
        p = subprocess.run([name] + list(args), stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=30)
        return (p.stdout or "").strip().splitlines()[0] if p.stdout else "?"
    except Exception:
        return "not found"


def _install_proxy(cfg, log):
    proxies = {}
    if cfg["network"].get("http_proxy"):
        proxies["http"] = cfg["network"]["http_proxy"]
    if cfg["network"].get("https_proxy"):
        proxies["https"] = cfg["network"]["https_proxy"]
    if proxies:
        opener = urllib.request.build_opener(urllib.request.ProxyHandler(proxies))
        urllib.request.install_opener(opener)
        log.info("Using HTTP proxy: %s", proxies)


def http_get_text(url: str, cfg, log) -> str:
    retries = int(cfg["network"]["retries"])
    delay = float(cfg["network"]["retry_delay_seconds"])
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "wheat-panproteome/2.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            last = e
            log.warning("GET text failed (%d/%d) %s: %s", attempt, retries, url, e)
            time.sleep(delay * attempt)
    raise RuntimeError(f"Failed to fetch {url}: {last}")


def http_download(url: str, dest: str, cfg, log, gunzip: bool = False) -> str:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        log.info("cache hit  : %s", os.path.basename(dest))
        return dest
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    retries = int(cfg["network"]["retries"])
    delay = float(cfg["network"]["retry_delay_seconds"])
    tmp = dest + (".gz.part" if gunzip and url.endswith(".gz") else ".part")
    last = None
    for attempt in range(1, retries + 1):
        try:
            log.info("download   : %s", url)
            req = urllib.request.Request(url, headers={"User-Agent": "wheat-panproteome/2.0"})
            with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
                shutil.copyfileobj(r, fh, length=1 << 20)
            if gunzip and url.endswith(".gz"):
                with gzip.open(tmp, "rb") as gi, open(dest, "wb") as go:
                    shutil.copyfileobj(gi, go, length=1 << 20)
                os.remove(tmp)
            else:
                os.replace(tmp, dest)
            if os.path.getsize(dest) == 0:
                raise RuntimeError("downloaded file is empty")
            return dest
        except Exception as e:
            last = e
            log.warning("download failed (%d/%d): %s", attempt, retries, e)
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except OSError:
                    pass
            time.sleep(delay * attempt)
    raise RuntimeError(f"Failed to download {url}: {last}")


def list_dir_hrefs(url: str, cfg, log):
    if not url.endswith("/"):
        url += "/"
    html = http_get_text(url, cfg, log)
    return re.findall(r'href="([^"?/][^"]*)"', html)


DEFAULT_CONFIG = {
    "whole_proteome": {
        "enabled": True,
        "resume": True,
        "work_dir": "",
        "collapse_isoforms": True,
        "min_identity": 40, "min_coverage": 60, "evalue": 1e-5,
        "cluster_novel": True, "cluster_approx_id": 90, "cluster_member_cover": 80,
        "sensitivity": "",
        "block_size_gb": 2.0,
        "index_chunks": 4,
        "memory_limit_gb": 48,
        "fallback_max_target_seqs": 25,
    },
    "project": {
        "name": "wheat_rca_cpn60_panproteome",
        "inputs_dir": "./inputs",
        "outputs_dir": "./outputs",
    },
    "network": {"http_proxy": "", "https_proxy": "", "retries": 4, "retry_delay_seconds": 5},
    "ensembl": {"release": "current",
                "division_ftp": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants",
                "division_ftp_fallbacks": [
                    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/plants",
                    "http://ftp.ensemblgenomes.org/pub/plants",
                    "https://ftp.ensemblgenomes.org/pub/plants",
                ]},
    "references": {
        "uniprot_rest": "https://rest.uniprot.org/uniprotkb/stream",
        "max_per_family": 200,
        "RCA": {
            "uniprot_query": '(taxonomy_id:4565) AND (protein_name:"rubisco activase" '
                             'OR protein_name:"ribulose bisphosphate carboxylase/oxygenase activase" '
                             'OR gene:rca)',
            "fasta": "",
        },
        "CPN60": {
            "uniprot_query": '(taxonomy_id:4565) AND (protein_name:"chaperonin 60" '
                             'OR protein_name:"60 kDa chaperonin" '
                             'OR protein_name:"chaperonin-60" '
                             'OR protein_name:"RuBisCO subunit binding-protein" '
                             'OR protein_name:"RuBisCO large subunit-binding protein" '
                             'OR protein_name:"Chaperonin CPN60" '
                             'OR gene:cpn60 OR gene:cpn60a OR gene:cpn60b OR gene:hsp60)',
            "fasta": "",
            "min_references": 6,
        },
    },
    "cultivars": {},
    "homology": {
        "min_identity": 40.0, "min_query_coverage": 60.0, "max_evalue": 1e-20,
        "diamond_sensitivity": "--more-sensitive", "max_target_seqs": 500, "threads": 8,
    },
    "alignment": {"mafft_args": "--auto --maxiterate 1000",
                  "collapse_support_below": 0.8,
                  "root": "midpoint",
                  "outgroup": {"RCA": [], "CPN60": []},
                  "add_outgroup_sequences": True},
    "thermostability": {"RCA_residues": [], "CPN60_residues": []},
    "report": {"draw_circos": False,
               "candidate_clade": []},
    "external": {
        "enabled": True,
        "annotate_proteotypes": True,
        "fetch_outgroups": True,
        "cache_dir": "",
        "max_accessions": 400,
        "outgroups": {
            "RCA": ["Q40071", "P10896", "Q9SB58"],
            "CPN60": ["P0A6F5"],
        },
    },
}


def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: str, log=None) -> dict:
    user = {}
    if path and os.path.exists(path):
        try:
            import yaml
        except Exception as e:
            raise RuntimeError("PyYAML is required to read config.yaml") from e
        with open(path) as fh:
            user = yaml.safe_load(fh) or {}
    elif path and log:
        log.warning("Config file %s not found; using built-in defaults.", path)
    cfg = _deep_merge(DEFAULT_CONFIG, user)
    if not cfg["cultivars"]:
        cfg["cultivars"] = {
            "chinese_spring": {"ensembl_strain": "triticum_aestivum", "group": "unknown"},
            "jagger":       {"ensembl_strain": "triticum_aestivum_jagger", "group": "unknown"},
            "julius":       {"ensembl_strain": "triticum_aestivum_julius", "group": "unknown"},
            "mace":         {"ensembl_strain": "triticum_aestivum_mace", "group": "unknown"},
            "arinalrfor":   {"ensembl_strain": "triticum_aestivum_arinalrfor", "group": "unknown"},
            "lancer":       {"ensembl_strain": "triticum_aestivum_lancer", "group": "unknown"},
            "landmark":     {"ensembl_strain": "triticum_aestivum_landmark", "group": "unknown"},
            "stanley":      {"ensembl_strain": "triticum_aestivum_stanley", "group": "unknown"},
            "norin61":      {"ensembl_strain": "triticum_aestivum_norin61", "group": "unknown"},
            "mattis":       {"ensembl_strain": "triticum_aestivum_mattis", "group": "unknown"},
        }
    return cfg


def ensembl_bases(cfg):
    raw = [cfg["ensembl"].get("division_ftp")] + list(
        cfg["ensembl"].get("division_ftp_fallbacks", []) or [])
    seen, out = set(), []
    for b in raw:
        if not b:
            continue
        b = b.rstrip("/")
        if b not in seen:
            seen.add(b)
            out.append(b)
    return out


def resolve_release(cfg, log):
    bases = ensembl_bases(cfg)
    rel = str(cfg["ensembl"]["release"]).strip()
    if rel and rel.lower() != "current":
        release_dir = rel if rel.startswith("release-") else f"release-{rel}"
        for base in bases:
            try:
                list_dir_hrefs(f"{base}/{release_dir}/fasta/", cfg, log)
                log.info("Using Ensembl host %s (%s)", base, release_dir)
                return base, release_dir
            except Exception as e:
                log.warning("Ensembl host not usable (%s): %s", base, e)
        log.warning("No host confirmed %s; proceeding with %s anyway.", release_dir, bases[0])
        return bases[0], release_dir
    for base in bases:
        try:
            hrefs = list_dir_hrefs(base, cfg, log)
            nums = [int(m.group(1)) for h in hrefs
                    for m in [re.match(r"release-(\d+)", h)] if m]
            if nums:
                release_dir = f"release-{max(nums)}"
                log.info("Resolved Ensembl host %s release %s", base, release_dir)
                return base, release_dir
        except Exception as e:
            log.warning("Could not list Ensembl host %s: %s", base, e)
    log.warning("Could not resolve a release from any host; trying 'current' on %s.", bases[0])
    return bases[0], "current"


def fetch_reference_sequences(cfg, refs_dir, log):
    os.makedirs(refs_dir, exist_ok=True)
    qid2family = {}
    per_family_paths = {}
    for fam in FAMILIES:
        spec = cfg["references"][fam]
        out = os.path.join(refs_dir, f"{fam}_ref.faa")
        if spec.get("fasta"):
            src = spec["fasta"]
            if not os.path.exists(src):
                raise RuntimeError(f"references.{fam}.fasta set but not found: {src}")
            shutil.copyfile(src, out)
            log.info("[ref] %s: using local FASTA %s", fam, src)
        else:
            q = spec["uniprot_query"]
            size = int(cfg["references"].get("max_per_family", 50))
            url = (cfg["references"]["uniprot_rest"]
                   + "?format=fasta&size=" + str(size)
                   + "&query=" + urllib.parse.quote(q))
            log.info("[ref] %s: querying UniProt", fam)
            http_download(url, out, cfg, log, gunzip=False)
        n = 0
        with open(out) as fh:
            for line in fh:
                if line.startswith(">"):
                    qid = line[1:].strip().split()[0]
                    qid2family[qid] = fam
                    n += 1
        if n == 0:
            raise RuntimeError(
                f"No reference sequences obtained for {fam}. Refine references.{fam}."
                f"uniprot_query or supply references.{fam}.fasta.")
        log.info("[ref] %s: %d reference sequence(s)", fam, n)
        per_family_paths[fam] = out

    combined = os.path.join(refs_dir, "references_combined.faa")
    with open(combined, "w") as out_fh:
        for fam in FAMILIES:
            with open(per_family_paths[fam]) as fh:
                shutil.copyfileobj(fh, out_fh)
    per_family_paths["combined"] = combined
    per_family_paths["qid2family"] = qid2family

    rows = []
    for fam in FAMILIES:
        n_fam = sum(1 for f in qid2family.values() if f == fam)
        min_ref = int(cfg["references"][fam].get("min_references", 0) or 0)
        if min_ref and n_fam < min_ref:
            log.warning("[ref] %s: only %d reference sequence(s) (expected >= %d). A small "
                        "query set recruits paralogues and superfamily members; broaden "
                        "references.%s.uniprot_query or supply a curated FASTA before "
                        "interpreting family statistics.", fam, n_fam, min_ref, fam)
        with open(per_family_paths[fam]) as fh:
            for line in fh:
                if line.startswith(">"):
                    parts = line[1:].strip().split(None, 1)
                    acc = parts[0].split("|")[1] if "|" in parts[0] else parts[0]
                    rows.append({"family": fam, "query_id": parts[0], "accession": acc,
                                 "description": parts[1] if len(parts) > 1 else ""})
    try:
        qtab = os.path.join(refs_dir, "reference_query_sets.csv")
        pd.DataFrame(rows).to_csv(qtab, index=False)
        per_family_paths["query_table"] = qtab
        log.info("[ref] query set table -> %s", os.path.basename(qtab))
    except Exception as e:
        log.warning("[ref] could not write reference_query_sets.csv: %s", e)
    return per_family_paths


def resolve_pep_url(base, release_dir, strain, cfg, log):
    pep_dir = f"{base}/{release_dir}/fasta/{strain}/pep/"
    hrefs = list_dir_hrefs(pep_dir, cfg, log)
    cands = [h for h in hrefs if re.search(r"\.pep\.all\.fa\.gz$", h)]
    if not cands:
        cands = [h for h in hrefs if re.search(r"\.pep\..*\.fa\.gz$", h)]
    if not cands:
        raise RuntimeError(f"No pep FASTA found under {pep_dir}")
    return pep_dir + sorted(cands)[0]


def fetch_cultivar_proteomes(cfg, inputs_dir, base, release_dir, log, only=None):
    prot_dir = os.path.join(inputs_dir, "proteomes")
    os.makedirs(prot_dir, exist_ok=True)
    result = {}
    for cult, meta in cfg["cultivars"].items():
        if only and cult not in only:
            continue
        local = meta.get("proteome_path", "")
        dest = os.path.join(prot_dir, f"{sanitize_name(cult)}.pep.fa")
        try:
            if local:
                if not os.path.exists(local):
                    raise RuntimeError(f"proteome_path not found: {local}")
                if local.endswith(".gz"):
                    if not (os.path.exists(dest) and os.path.getsize(dest) > 0):
                        with gzip.open(local, "rb") as gi, open(dest, "wb") as go:
                            shutil.copyfileobj(gi, go)
                else:
                    dest = local
                log.info("[proteome] %-14s local -> %s", cult, dest)
            else:
                strain = meta.get("ensembl_strain")
                if not strain:
                    raise RuntimeError("neither proteome_path nor ensembl_strain set")
                url = resolve_pep_url(base, release_dir, strain, cfg, log)
                http_download(url, dest, cfg, log, gunzip=True)
                log.info("[proteome] %-14s %s", cult, os.path.basename(dest))
            if os.path.getsize(dest) == 0:
                raise RuntimeError("empty proteome file")
            result[cult] = dest
        except Exception as e:
            log.warning("[proteome] %-14s SKIPPED: %s", cult, e)
    if not result:
        raise RuntimeError("No cultivar proteomes could be obtained. Check network/proxy, "
                           "the Ensembl release, or supply proteome_path entries in config.")
    log.info("Obtained %d/%d cultivar proteomes.", len(result), len(cfg["cultivars"]))
    return result


DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                  "qcovhsp", "scovhsp"]


def diamond_search(diamond, references, proteome_fa, work_dir, cfg, log):
    os.makedirs(work_dir, exist_ok=True)
    threads = str(cfg["homology"]["threads"])
    db = os.path.join(work_dir, "proteome_db")
    run_cmd([diamond, "makedb", "--in", proteome_fa, "--db", db,
             "--threads", threads, "--quiet"], log)
    out = os.path.join(work_dir, "hits.tsv")
    cmd = [diamond, "blastp", "--query", references, "--db", db + ".dmnd",
           "--out", out, "--outfmt", "6", *DIAMOND_FIELDS,
           "--evalue", str(cfg["homology"]["max_evalue"]),
           "--id", str(cfg["homology"]["min_identity"]),
           "--query-cover", str(cfg["homology"]["min_query_coverage"]),
           "--max-target-seqs", str(cfg["homology"]["max_target_seqs"]),
           "--threads", threads, "--quiet"]
    sens = cfg["homology"].get("diamond_sensitivity", "").strip()
    if sens:
        cmd.append(sens)
    run_cmd(cmd, log)
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        return pd.DataFrame(columns=DIAMOND_FIELDS)
    df = pd.read_csv(out, sep="\t", names=DIAMOND_FIELDS)
    df = df[(df["pident"] >= cfg["homology"]["min_identity"])
            & (df["qcovhsp"] >= cfg["homology"]["min_query_coverage"])
            & (df["evalue"] <= cfg["homology"]["max_evalue"])]
    return df.reset_index(drop=True)


def assign_families(hits: pd.DataFrame, qid2family: dict):
    if hits.empty:
        return {}
    hits = hits.copy()
    hits["family"] = hits["qseqid"].map(lambda q: qid2family.get(q, "NA"))
    hits = hits[hits["family"] != "NA"]
    best = {}
    for sid, grp in hits.sort_values(["evalue", "bitscore"],
                                     ascending=[True, False]).groupby("sseqid"):
        top = grp.iloc[0]
        best[sid] = (top["family"], float(top["pident"]), float(top["evalue"]))
    return best


def extract_members(cultivar, proteome_fa, subject_family, log):
    if not _HAVE_BIO:
        raise RuntimeError("Biopython required for sequence extraction.")
    members = []
    if not subject_family:
        return members
    idx = SeqIO.index(proteome_fa, "fasta")
    try:
        for sid, (fam, pid, ev) in subject_family.items():
            if sid not in idx:
                alt = next((k for k in idx.keys() if k.split(".")[0] == sid.split(".")[0]), None)
                if alt is None:
                    log.debug("subject %s not found in %s index", sid, cultivar)
                    continue
                sid = alt
            seq = str(idx[sid].seq).rstrip("*")
            if not seq:
                continue
            gene_id = re.sub(r"\.\d+$", "", sid)
            _mc = re.match(r"^Traes[A-Z0-9]{2,4}(\d[ABD]|U)", sid)
            chrom = "un" if (not _mc or _mc.group(1) == "U") else _mc.group(1)
            members.append({"cultivar": cultivar, "family": fam, "source_id": sid,
                            "gene": gene_id, "chromosome": chrom,
                            "subgenome": ("un" if chrom == "un" else chrom[-1]),
                            "homoeologous_group": ("un" if chrom == "un" else chrom[0]),
                            "pident_to_ref": round(pid, 2), "evalue_to_ref": ev,
                            "length": len(seq), "seq": seq})
    finally:
        idx.close()
    return members


def dedupe_to_proteotypes(members, reference_cultivar="chinese_spring"):
    proteotypes = []
    member_records = []
    for fam in FAMILIES:
        fam_members = [m for m in members if m["family"] == fam]
        by_seq = defaultdict(list)
        for m in fam_members:
            by_seq[m["seq"]].append(m)
        def sort_key(item):
            seq, ms = item
            culs = {x["cultivar"] for x in ms}
            return (-len(culs), -len(seq), seq)
        for i, (seq, ms) in enumerate(sorted(by_seq.items(), key=sort_key), start=1):
            pt_id = f"{fam}_PT{i:03d}"
            culs = sorted({x["cultivar"] for x in ms})
            in_ref = reference_cultivar in culs
            proteotypes.append({
                "proteotype_id": pt_id, "family": fam, "length": len(seq),
                "n_cultivars": len(culs), "cultivars": ";".join(culs),
                "present_in_reference": in_ref,
                "novel_vs_reference": (not in_ref),
                "source_ids": ";".join(sorted(x["source_id"] for x in ms)),
                "seq": seq,
            })
            for m in ms:
                rec = dict(m)
                rec["proteotype_id"] = pt_id
                member_records.append(rec)
    return proteotypes, member_records


def build_pav_matrix(proteotypes, cultivars):
    rows = []
    for pt in proteotypes:
        present = set(pt["cultivars"].split(";")) if pt["cultivars"] else set()
        row = {"proteotype_id": pt["proteotype_id"], "family": pt["family"]}
        for c in cultivars:
            row[c] = int(c in present)
        rows.append(row)
    return pd.DataFrame(rows)


def write_family_fasta(proteotypes, fam, path):
    if not _HAVE_BIO:
        raise RuntimeError("Biopython required to write FASTA.")
    recs = [SeqRecord(Seq(pt["seq"]), id=pt["proteotype_id"],
                      description=f"{fam} novel={pt['novel_vs_reference']} "
                                  f"n_cultivars={pt['n_cultivars']}")
            for pt in proteotypes if pt["family"] == fam]
    with open(path, "w") as fh:
        SeqIO.write(recs, fh, "fasta")
    return len(recs)


def run_mafft(mafft, in_fa, out_aln, cfg, log):
    args = cfg["alignment"]["mafft_args"].split()
    with open(out_aln, "w") as fh:
        proc = subprocess.run([mafft, *args, in_fa], stdout=fh,
                              stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        log.error("MAFFT failed: %s", (proc.stderr or "")[-500:])
        raise RuntimeError("MAFFT alignment failed")
    return out_aln


def run_fasttree(fasttree, in_aln, out_nwk, log):
    with open(out_nwk, "w") as fh:
        proc = subprocess.run([fasttree, "-quiet", in_aln], stdout=fh,
                              stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        log.error("FastTree failed: %s", (proc.stderr or "")[-500:])
        raise RuntimeError("FastTree failed")
    return out_nwk


def plot_circular_tree(nwk, out_png, title, log, novel_ids=None, family=""):
    novel_ids = novel_ids or set()
    try:
        tree = Phylo.read(nwk, "newick")
        tree.ladderize()
        terminals = tree.get_terminals()
        n = len(terminals)
        big = n > 60
        if n < 3:
            log.warning("Tree for %s has <3 tips; skipping circular tree.", title)
            return None

        gap = math.radians(14)
        span = 2 * math.pi - gap
        start = math.pi / 2 + gap / 2
        leaf_angle = {t: start + span * i / (n - 1) for i, t in enumerate(terminals)}

        depths = tree.depths()
        if not depths or max(depths.values()) == 0:
            depths = tree.depths(unit_branch_lengths=True)
        maxd = max(depths.values()) or 1.0
        radius = {c: depths[c] / maxd for c in depths}

        angle = {}
        def _angle(cl):
            if cl.is_terminal():
                angle[cl] = leaf_angle[cl]
            else:
                kids = [_angle(c) for c in cl.clades]
                angle[cl] = sum(kids) / len(kids)
            return angle[cl]
        _angle(tree.root)

        def xy(r, a):
            return r * math.cos(a), r * math.sin(a)

        fig, ax = plt.subplots(figsize=((14, 14) if big else (9.5, 9.5)))
        ax.set_aspect("equal"); ax.axis("off")
        for parent in tree.find_clades(order="level"):
            if not parent.clades:
                continue
            rp = radius[parent]
            child_angles = [angle[c] for c in parent.clades]
            arc = np.linspace(min(child_angles), max(child_angles), 60)
            ax.plot(rp * np.cos(arc), rp * np.sin(arc), color="#444", lw=0.6)
            for c in parent.clades:
                x0, y0 = xy(rp, angle[c]); x1, y1 = xy(radius[c], angle[c])
                ax.plot([x0, x1], [y0, y1], color="#444", lw=0.6)

        r_label = 1.06
        for t in terminals:
            a = leaf_angle[t]; r = radius[t]
            tx, ty = xy(r, a)
            col = FIG_RED if t.name in novel_ids else FIG_BLUE
            ax.plot([tx], [ty], marker="o", ms=2.6, color=col, zorder=3)
            gx, gy = xy(r_label, a)
            ax.plot([tx, gx], [ty, gy], color="#dddddd", lw=0.4, zorder=1)
            deg = math.degrees(a) % 360
            lab = t.name.split("_", 1)[1] if "_" in t.name else t.name
            if 90 < deg < 270:
                rot, ha = deg + 180, "right"
            else:
                rot, ha = deg, "left"
            lx, ly = xy(r_label + 0.015, a)
            ax.text(lx, ly, lab, rotation=rot, rotation_mode="anchor",
                    ha=ha, va="center", fontsize=(3.2 if big else 6.3))

        _legend_below(fig, ax,
                      [Patch(facecolor=FIG_BLUE, label="present in Chinese Spring"),
                       Patch(facecolor=FIG_RED, label="novel vs Chinese Spring")],
                      ncol=2, fontsize=9)
        _lim = 1.62 if big else 1.5
        ax.set_xlim(-_lim, _lim); ax.set_ylim(-_lim, _lim)
        fig.savefig(out_png, dpi=600, bbox_inches="tight")
        plt.close(fig)
        return out_png
    except Exception as e:
        log.warning("Could not render circular tree for %s: %s", title, e)
        return None


def sanitize_aa(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in STANDARD_AA)


def aliphatic_index(seq: str) -> float:
    n = len(seq)
    if n == 0:
        return float("nan")
    c = Counter(seq)
    a = c["A"] / n * 100.0
    v = c["V"] / n * 100.0
    il = (c["I"] + c["L"]) / n * 100.0
    return a + 2.9 * v + 3.9 * il


def gravy(seq: str) -> float:
    if not seq:
        return float("nan")
    return sum(KD[c] for c in seq) / len(seq)


def molecular_weight(seq: str) -> float:
    if not seq:
        return float("nan")
    return sum(AA_MASS[c] for c in seq) + WATER


def aromaticity(seq: str) -> float:
    if not seq:
        return float("nan")
    return sum(seq.count(x) for x in "FWY") / len(seq)


def _charge_at_pH(counts, pH):
    pos = 1.0 / (1 + 10 ** (pH - PK["Nterm"]))
    for aa in ("K", "R", "H"):
        pos += counts.get(aa, 0) / (1 + 10 ** (pH - PK[aa]))
    neg = 1.0 / (1 + 10 ** (PK["Cterm"] - pH))
    for aa in ("D", "E", "C", "Y"):
        neg += counts.get(aa, 0) / (1 + 10 ** (PK[aa] - pH))
    return pos - neg


def isoelectric_point(seq: str) -> float:
    if not seq:
        return float("nan")
    counts = Counter(seq)
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        ch = _charge_at_pH(counts, mid)
        if abs(ch) < 1e-4:
            return round(mid, 2)
        if ch > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 2)


def physicochemical_features(seq_raw: str) -> dict:
    seq = sanitize_aa(seq_raw)
    n = len(seq)
    feats = {"length": len(seq_raw), "length_std_aa": n}
    if n == 0:
        return feats
    pct = {aa: seq.count(aa) / n for aa in STANDARD_AA}
    feats.update({
        "mw": round(molecular_weight(seq), 2),
        "pI": isoelectric_point(seq),
        "aromaticity": round(aromaticity(seq), 4),
        "gravy": round(gravy(seq), 4),
        "aliphatic_index": round(aliphatic_index(seq), 2),
        "frac_acidic_DE": round(pct["D"] + pct["E"], 4),
        "frac_basic_KR": round(pct["K"] + pct["R"], 4),
        "frac_charged_DEKRH": round(sum(pct[a] for a in "DEKRH"), 4),
        "net_charge_frac": round((pct["K"] + pct["R"]) - (pct["D"] + pct["E"]), 4),
        "frac_hydrophobic": round(sum(pct[a] for a in HYDROPHOBIC), 4),
        "frac_proline": round(pct["P"], 4),
        "frac_glycine": round(pct["G"], 4),
        "frac_arg_of_RK": round(pct["R"] / (pct["R"] + pct["K"]), 4) if (pct["R"] + pct["K"]) else float("nan"),
    })
    if _HAVE_PROTPARAM:
        try:
            pa = ProteinAnalysis(seq)
            feats["instability_index"] = round(pa.instability_index(), 2)
            helix, turn, sheet = pa.secondary_structure_fraction()
            feats["frac_helix"] = round(helix, 4)
            feats["frac_sheet"] = round(sheet, 4)
        except Exception:
            feats["instability_index"] = float("nan")
    else:
        feats["instability_index"] = float("nan")
    return feats


THERMO_FEATURES = ["aliphatic_index", "gravy", "instability_index", "frac_charged_DEKRH",
                   "frac_acidic_DE", "frac_basic_KR", "net_charge_frac", "frac_proline",
                   "frac_glycine", "frac_arg_of_RK", "pI", "aromaticity"]

GROUP_COMPARISON_COLS = ["family", "feature", "n_resistant", "n_sensitive",
                         "median_resistant", "median_sensitive", "U", "p_value", "q_value_BH"]
KNOWN_RESIDUE_COLS = ["family", "reference", "ref_position", "aln_column", "proteotype_id",
                      "residue", "cultivars", "note"]
POLY_COLS = ["family", "aln_column", "n_variants", "residues", "entropy_bits", "gap_fraction"]
POS_ASSOC_COLS = ["family", "aln_column", "p_value", "resistant_dist", "sensitive_dist",
                  "q_value_BH"]


def benjamini_hochberg(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def groups_from_config(cfg):
    resistant, sensitive = [], []
    for c, meta in cfg["cultivars"].items():
        g = str(meta.get("group", "unknown")).lower()
        if g == "resistant":
            resistant.append(c)
        elif g == "sensitive":
            sensitive.append(c)
    return resistant, sensitive


def compare_group_features(member_features_df, resistant, sensitive, log):
    if not _HAVE_SCIPY:
        log.warning("scipy unavailable; skipping group statistics.")
        return pd.DataFrame(columns=GROUP_COMPARISON_COLS)
    df = member_features_df.copy()
    df["group"] = df["cultivar"].map(
        lambda c: "resistant" if c in resistant else ("sensitive" if c in sensitive else "other"))
    df = df[df["group"].isin(["resistant", "sensitive"])]
    rows = []
    for fam in FAMILIES:
        sub = df[df["family"] == fam]
        for feat in THERMO_FEATURES:
            if feat not in sub.columns:
                continue
            r = sub.loc[sub["group"] == "resistant", feat].dropna().values
            s = sub.loc[sub["group"] == "sensitive", feat].dropna().values
            if len(r) < 2 or len(s) < 2:
                continue
            try:
                u, p = scipy_stats.mannwhitneyu(r, s, alternative="two-sided")
            except ValueError:
                u, p = float("nan"), float("nan")
            rows.append({"family": fam, "feature": feat,
                         "n_resistant": len(r), "n_sensitive": len(s),
                         "median_resistant": round(float(np.median(r)), 4),
                         "median_sensitive": round(float(np.median(s)), 4),
                         "U": u, "p_value": p})
    out = pd.DataFrame(rows, columns=(None if rows else GROUP_COMPARISON_COLS))
    if not out.empty:
        out["q_value_BH"] = benjamini_hochberg(out["p_value"].fillna(1.0).values)
    return out


def polymorphic_sites(aln_records, fam):
    ids = list(aln_records)
    if len(ids) < 2:
        return pd.DataFrame(columns=["family", "aln_column", "n_variants",
                                     "residues", "entropy_bits", "gap_fraction"])
    L = len(aln_records[ids[0]])
    rows = []
    for col in range(L):
        column = [aln_records[i][col] for i in ids]
        counts = Counter(column)
        nongap = {r: c for r, c in counts.items() if r != "-"}
        if len(set(column)) <= 1:
            continue
        total = len(column)
        ent = -sum((c / total) * math.log2(c / total) for c in counts.values())
        gapf = counts.get("-", 0) / total
        dist = ",".join(f"{r}:{c}" for r, c in sorted(counts.items()))
        rows.append({"family": fam, "aln_column": col + 1,
                     "n_variants": len(nongap), "residues": dist,
                     "entropy_bits": round(ent, 3), "gap_fraction": round(gapf, 3)})
    return pd.DataFrame(rows, columns=(None if rows else POLY_COLS))


def position_group_association(aln_records, member_records, resistant, sensitive, fam, log):
    if not _HAVE_SCIPY:
        return pd.DataFrame(columns=POS_ASSOC_COLS)
    ids = list(aln_records)
    if len(ids) < 2:
        return pd.DataFrame(columns=POS_ASSOC_COLS)
    L = len(aln_records[ids[0]])
    labelled = []
    for m in member_records:
        if m["family"] != fam:
            continue
        g = ("resistant" if m["cultivar"] in resistant
             else "sensitive" if m["cultivar"] in sensitive else None)
        if g is None:
            continue
        aseq = aln_records.get(m["proteotype_id"])
        if aseq is None:
            continue
        labelled.append((g, aseq))
    if sum(1 for g, _ in labelled if g == "resistant") < 2 or \
       sum(1 for g, _ in labelled if g == "sensitive") < 2:
        return pd.DataFrame(columns=POS_ASSOC_COLS)
    rows = []
    for col in range(L):
        residues = [(g, aseq[col]) for g, aseq in labelled]
        symbols = sorted({r for _, r in residues})
        if len(symbols) <= 1:
            continue
        table = np.zeros((2, len(symbols)), dtype=int)
        gidx = {"resistant": 0, "sensitive": 1}
        sidx = {s: i for i, s in enumerate(symbols)}
        for g, r in residues:
            table[gidx[g], sidx[r]] += 1
        table = table[:, table.sum(axis=0) > 0]
        if table.shape[1] < 2 or (table.sum(axis=1) < 2).any():
            continue
        try:
            _, p, _, _ = scipy_stats.chi2_contingency(table)
        except Exception:
            continue
        dist_r = ",".join(f"{s}:{table[0, i]}" for i, s in enumerate(symbols))
        dist_s = ",".join(f"{s}:{table[1, i]}" for i, s in enumerate(symbols))
        rows.append({"family": fam, "aln_column": col + 1, "p_value": p,
                     "resistant_dist": dist_r, "sensitive_dist": dist_s})
    out = pd.DataFrame(rows, columns=(None if rows else POS_ASSOC_COLS))
    if not out.empty:
        out["q_value_BH"] = benjamini_hochberg(out["p_value"].values)
        out = out.sort_values("p_value").reset_index(drop=True)
    return out


def map_known_residues(cfg, aln_by_family, member_records, log):
    rows = []
    for fam in FAMILIES:
        entries = cfg["thermostability"].get(f"{fam}_residues", []) or []
        if not entries:
            continue
        aln = aln_by_family.get(fam, {})
        if not aln:
            continue
        pt2culs = defaultdict(set)
        for m in member_records:
            if m["family"] == fam:
                pt2culs[m["proteotype_id"]].add(m["cultivar"])
        for entry in entries:
            ref_id = entry.get("ref")
            pos = int(entry.get("pos", 0))
            note = entry.get("note", "")
            if ref_id not in aln:
                log.warning("[known-residue] %s ref '%s' not among aligned proteotypes; skipped.",
                            fam, ref_id)
                continue
            aseq = aln[ref_id]
            col = None
            seen = 0
            for i, ch in enumerate(aseq):
                if ch != "-":
                    seen += 1
                    if seen == pos:
                        col = i
                        break
            if col is None:
                log.warning("[known-residue] %s pos %d out of range for %s.", fam, pos, ref_id)
                continue
            for pt, aseq2 in aln.items():
                rows.append({"family": fam, "reference": ref_id, "ref_position": pos,
                             "aln_column": col + 1, "proteotype_id": pt,
                             "residue": aseq2[col],
                             "cultivars": ";".join(sorted(pt2culs.get(pt, []))),
                             "note": note})
    return pd.DataFrame(rows, columns=(None if rows else KNOWN_RESIDUE_COLS))


def write_sqlite(db_path, tables, log):
    conn = sqlite3.connect(db_path)
    try:
        for name, df in tables.items():
            if df is None:
                continue
            if getattr(df, "shape", (0, 0))[1] == 0:
                log.warning("SQLite: skipping table %s (no columns).", name)
                continue
            df.to_sql(name, conn, if_exists="replace", index=False)
        conn.commit()
    finally:
        conn.close()
    log.info("Wrote SQLite catalog: %s", db_path)


def plot_pav_heatmap(pav_df, fam, out_png, log):
    sub = pav_df[pav_df["family"] == fam]
    if sub.empty:
        return None
    mat = sub.drop(columns=["family"]).set_index("proteotype_id")
    if mat.shape[0] == 0 or mat.shape[1] == 0:
        return None
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * mat.shape[1] + 3),
                                    min(30.0, max(3, 0.4 * mat.shape[0] + 2))))
    if _HAVE_SNS:
        sns.heatmap(mat, cmap=["#f0f0f0", "#2166ac"], cbar=False, linewidths=0.5,
                    linecolor="white", ax=ax, vmin=0, vmax=1)
    else:
        ax.imshow(mat.values, aspect="auto", cmap="Blues")
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(mat.columns, rotation=90)
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index)
    ax.set_xlabel("cultivar")
    ax.set_ylabel("proteotype")
    fig.tight_layout()
    fig.savefig(out_png, dpi=600)
    plt.close(fig)
    return out_png


def plot_feature_boxplots(member_features_df, resistant, sensitive, fam, out_png, log):
    df = member_features_df[member_features_df["family"] == fam].copy()
    df["group"] = df["cultivar"].map(
        lambda c: "resistant" if c in resistant else ("sensitive" if c in sensitive else "other"))
    df = df[df["group"].isin(["resistant", "sensitive"])]
    feats = [f for f in ["aliphatic_index", "gravy", "instability_index",
                         "frac_charged_DEKRH", "net_charge_frac", "pI"]
             if f in df.columns and df[f].notna().any()]
    if df.empty or not feats:
        return None
    ncol = 3
    nrow = int(math.ceil(len(feats) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.2 * nrow), squeeze=False)
    for k, feat in enumerate(feats):
        ax = axes[k // ncol][k % ncol]
        data = [df.loc[df["group"] == g, feat].dropna().values
                for g in ("resistant", "sensitive")]
        try:
            ax.boxplot(data, tick_labels=["resistant", "sensitive"], showmeans=True)
        except TypeError:
            ax.boxplot(data, labels=["resistant", "sensitive"], showmeans=True)
        ax.set_title(feat)
    for k in range(len(feats), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.tight_layout()
    fig.savefig(out_png, dpi=600)
    plt.close(fig)
    return out_png


def plot_entropy(poly_df, fam, out_png, log):
    sub = poly_df[poly_df["family"] == fam].sort_values("aln_column")
    if sub.empty:
        return None
    disp = FAM_DISPLAY.get(fam, fam)
    x = sub["aln_column"].to_numpy()
    y = sub["entropy_bits"].to_numpy()
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.bar(x, y, width=1.0, color=FIG_RED, alpha=0.85, label="per-column entropy")
    if len(y) >= 15:
        w = max(5, len(y) // 40)
        trend = pd.Series(y).rolling(w, center=True, min_periods=1).mean()
        ax.plot(x, trend, color="#1a1a1a", lw=1.6, label=f"rolling mean (w={w})")
    ax.set_xlabel("Alignment column", fontsize=14)
    ax.set_ylabel("Shannon entropy (bits)", fontsize=14)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_ylim(0, max(1.0, y.max() * 1.32))
    ax.margins(x=0.01)
    ax.spines[["top", "right"]].set_visible(False)
    imax = int(np.argmax(y))
    frac = (x[imax] - x.min()) / max(1e-9, (x.max() - x.min()))
    ha = "right" if frac > 0.75 else ("left" if frac < 0.25 else "center")
    ax.annotate(f"max {y[imax]:.2f} bits (col {int(x[imax])})",
                xy=(x[imax], y[imax]), xytext=(0, 14), textcoords="offset points",
                ha=ha, va="bottom", fontsize=10,
                arrowprops=dict(arrowstyle="->", lw=1.0, color="#1a1a1a"))
    ax.legend(fontsize=10, frameon=False, loc="lower left",
              bbox_to_anchor=(0.0, 1.005), ncol=2, borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _novel_set(pav_df, fam):
    sub = pav_df[pav_df["family"] == fam]
    if sub.empty or "chinese_spring" not in sub.columns:
        return set()
    return set(sub.loc[sub["chinese_spring"] == 0, "proteotype_id"])


def plot_circos_pav(pav_df, fam, out_png, log):
    sub = pav_df[pav_df["family"] == fam]
    if sub.empty:
        return None
    mat = sub.drop(columns=["family"]).set_index("proteotype_id")
    proteotypes = list(mat.index)
    cultivars = list(mat.columns)
    if not proteotypes or not cultivars:
        return None
    novel = _novel_set(pav_df, fam)
    disp = FAM_DISPLAY.get(fam, fam)

    n = len(proteotypes)
    gap = math.radians(26)
    span = 2 * math.pi - gap
    start = math.pi / 2 + gap / 2
    theta = start + (np.arange(n) + 0.5) * (span / n)
    width = (span / n) * 0.9
    r0, ring_h = 2.5, 1.0

    size = 11 if n <= 60 else 13
    fig, ax = plt.subplots(figsize=(size, size), subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("E")
    for i, cult in enumerate(cultivars):
        bottom = r0 + i * ring_h
        colors = []
        for pt in proteotypes:
            if mat.loc[pt, cult] == 1:
                colors.append(FIG_RED if pt in novel else FIG_BLUE)
            else:
                colors.append(FIG_ABSENT)
        ax.bar(theta, [ring_h * 0.9] * n, width=width, bottom=bottom,
               color=colors, edgecolor="white", linewidth=0.3, alpha=0.97)
        ax.text(math.pi / 2, bottom + ring_h * 0.45, clabel_display(cult),
                ha="center", va="center", fontsize=8)

    r_label = r0 + len(cultivars) * ring_h + 0.35
    for th, pt in zip(theta, proteotypes):
        deg = math.degrees(th) % 360
        lab = pt.split("_", 1)[1] if "_" in pt else pt
        if 90 < deg < 270:
            rot, ha = deg + 180, "right"
        else:
            rot, ha = deg, "left"
        ax.text(th, r_label, lab, rotation=rot, rotation_mode="anchor",
                ha=ha, va="center", fontsize=5 if n > 60 else 6.5)

    ax.set_xticks([]); ax.set_yticks([])
    ax.set_ylim(0, r_label + 0.9)
    ax.grid(False)
    ax.spines["polar"].set_visible(False)
    ax.legend(handles=[Patch(facecolor=FIG_BLUE, label="present (in Chinese Spring)"),
                       Patch(facecolor=FIG_RED, label="present, novel vs Chinese Spring"),
                       Patch(facecolor=FIG_ABSENT, edgecolor="#cccccc", label="absent")],
              loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=3,
              frameon=False, fontsize=10)
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _legend_below(fig, axes, handles, ncol=2, fontsize=10, pad=0.014):
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


def clabel_display(c):
    special = {"chinese_spring": "Chinese Spring", "arinalrfor": "ArinaLrFor",
               "norin61": "Norin61", "mattis": "SY Mattis", "stanley": "CDC Stanley"}
    return special.get(c, c.replace("_", " ").title())


def plot_overview(proteotypes, pav_df, cultivars, out_png, log):
    fams = list(FAMILIES)
    total = {f: sum(1 for p in proteotypes if p["family"] == f) for f in fams}
    noveln = {f: sum(1 for p in proteotypes if p["family"] == f and p["novel_vs_reference"])
              for f in fams}
    novel_ids = {f: _novel_set(pav_df, f) for f in fams}

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.5, 3.6))
    x = np.arange(len(fams)); w = 0.38
    axA.bar(x - w/2, [total[f] for f in fams], w, label="Total proteotypes",
            color=FIG_BLUE, edgecolor="white", linewidth=0.5)
    axA.bar(x + w/2, [noveln[f] for f in fams], w, label="Novel vs Chinese Spring",
            color=FIG_RED, edgecolor="white", linewidth=0.5)
    for xi, f in enumerate(fams):
        axA.text(xi - w/2, total[f] + 2, str(total[f]), ha="center", va="bottom", fontsize=8)
        axA.text(xi + w/2, noveln[f] + 2,
                 f"{noveln[f]} ({100*noveln[f]/max(1,total[f]):.0f}%)",
                 ha="center", va="bottom", fontsize=7.5)
    axA.set_xticks(x); axA.set_xticklabels([FAM_DISPLAY[f] for f in fams])
    axA.set_ylabel("Non-redundant proteotypes")
    axA.set_ylim(0, max(total.values()) * 1.25)
    axA.legend(frameon=False, fontsize=8, loc="upper left")
    axA.spines[["top", "right"]].set_visible(False)
    axA.set_title("a", loc="left", fontweight="bold", fontsize=12)

    counts = {f: {} for f in fams}
    for f in fams:
        pv = pav_df[(pav_df["family"] == f) & (pav_df["proteotype_id"].isin(novel_ids[f]))]
        for c in cultivars:
            counts[f][c] = int(pv[c].sum()) if not pv.empty else 0
    order = sorted(cultivars, key=lambda c: -(counts["RCA"][c] + counts["CPN60"][c]))
    xb = np.arange(len(order)); wb = 0.4
    axB.bar(xb - wb/2, [counts["RCA"][c] for c in order], wb, label="RCA",
            color=FIG_BLUE, edgecolor="white", linewidth=0.5)
    axB.bar(xb + wb/2, [counts["CPN60"][c] for c in order], wb, label="Cpn60",
            color="#E0A458", edgecolor="white", linewidth=0.5)
    axB.set_xticks(xb)
    axB.set_xticklabels([clabel_display(c) for c in order], rotation=45, ha="right", fontsize=7.5)
    axB.set_ylabel("Novel proteotypes carried")
    axB.legend(frameon=False, fontsize=8, title=None)
    axB.spines[["top", "right"]].set_visible(False)
    axB.set_title("b", loc="left", fontweight="bold", fontsize=12)

    fig.tight_layout()
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out_png


def plot_novel_pav(proteotypes, pav_df, fam, cultivars, out_png, log, max_rows=None):
    novel_ids = _novel_set(pav_df, fam)
    if not novel_ids:
        return None
    sub = pav_df[pav_df["proteotype_id"].isin(novel_ids)].copy()
    sub["prevalence"] = sub[cultivars].sum(axis=1)
    sub = sub.sort_values(["prevalence", "proteotype_id"], ascending=[False, True])
    n_total = len(sub)
    if max_rows and n_total > max_rows:
        sub = sub.head(max_rows)
    ids = sub["proteotype_id"].tolist()
    M = sub[cultivars].to_numpy(dtype=int)
    disp = FAM_DISPLAY.get(fam, fam)

    fig, ax = plt.subplots(figsize=(max(5.5, 0.55 * len(cultivars) + 2),
                                    max(3.5, 0.22 * len(ids) + 1.5)))
    ax.imshow(M, cmap=ListedColormap([FIG_ABSENT, FIG_BLUE]), aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(cultivars)))
    ax.set_xticklabels([clabel_display(c) for c in cultivars], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(ids)))
    ax.set_yticklabels([i.replace("_PT", " ") for i in ids], fontsize=6.5)
    ax.set_xticks(np.arange(-.5, len(cultivars), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(ids), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="both", length=0)
    shown = f"all {n_total}" if (not max_rows or n_total <= max_rows) else f"{len(ids)} of {n_total}"
    ax.set_ylabel(f"{disp} novel proteotypes ({shown})", fontsize=9)
    _legend_below(fig, ax,
                  [Patch(facecolor=FIG_BLUE, label="present"),
                   Patch(facecolor=FIG_ABSENT, edgecolor="#cccccc", label="absent")],
                  ncol=2, fontsize=8)
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _fmt_pct(n, d):
    return f"{(100.0 * n / d):.1f}%" if d else "n/a"


def _pangenome_partition(pav_df, cult_cols):
    M = pav_df[cult_cols].to_numpy(bool); prev = M.sum(1); N = len(cult_cols)
    soft = math.ceil(0.9 * N)
    return ({"core": int((prev == N).sum()),
             "soft-core": int(((prev >= soft) & (prev < N)).sum()),
             "shell": int(((prev >= 2) & (prev < soft)).sum()),
             "cloud": int((prev == 1).sum())}, prev)


def _rarefaction(pav_df, cult_cols, nperm=50, seed=0):
    M = pav_df[cult_cols].to_numpy(bool); N = len(cult_cols)
    rng = np.random.default_rng(seed)
    pan = np.zeros((nperm, N)); core = np.zeros((nperm, N)); idx = np.arange(N)
    for pp in range(nperm):
        rng.shuffle(idx); ca = np.zeros(M.shape[0], bool); cl = np.ones(M.shape[0], bool)
        for k in range(N):
            c = M[:, idx[k]]; ca |= c; cl &= c; pan[pp, k] = ca.sum(); core[pp, k] = cl.sum()
    return pan, core


def plot_pangenome_overview(pav_df, cult_cols, out_png, log,
                            title="Whole-proteome pan-genome (Chinese Spring reference)"):
    try:
        part, _ = _pangenome_partition(pav_df, cult_cols)
        pan, core = _rarefaction(pav_df, cult_cols); N = len(cult_cols)
        fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 3.9))
        cats = ["core", "soft-core", "shell", "cloud"]; vals = [part[c] for c in cats]
        tot = max(1, sum(vals))
        axA.bar(range(4), vals, color=[FIG_BLUE, "#6C9BD1", FIG_AMBER, FIG_RED], edgecolor="white")
        for i, v in enumerate(vals):
            axA.text(i, v + tot * 0.012, f"{v}\n({100 * v / tot:.0f}%)", ha="center", va="bottom", fontsize=8)
        axA.set_xticks(range(4))
        axA.set_xticklabels([f"Core\n(all {N})", "Soft-core\n(\u226590%)", "Shell", "Cloud\n(private)"], fontsize=8.5)
        axA.set_ylabel("Proteotype clusters", fontsize=10); axA.set_ylim(0, max(vals) * 1.22)
        axA.spines[["top", "right"]].set_visible(False)
        _panel_letter(axA, "a")
        ks = np.arange(1, N + 1)
        for arr, lab, c in [(pan, "Pan-proteome", FIG_BLUE), (core, "Core proteome", FIG_RED)]:
            axB.plot(ks, arr.mean(0), "-o", color=c, ms=3.5, label=lab)
            axB.fill_between(ks, arr.min(0), arr.max(0), color=c, alpha=0.15)
        axB.set_xlabel("Cultivars sampled", fontsize=10); axB.set_ylabel("Proteotype clusters", fontsize=10)
        axB.set_xticks(ks); axB.legend(frameon=False, fontsize=9)
        axB.spines[["top", "right"]].set_visible(False)
        _panel_letter(axB, "b")
        fig.tight_layout()
        fig.savefig(out_png, bbox_inches="tight", dpi=600); plt.close(fig)
        return out_png
    except Exception as e:
        log.warning("pangenome overview failed: %s", e); return None


def plot_prevalence_hist(pav_df, cult_cols, out_png, log):
    try:
        _, prev = _pangenome_partition(pav_df, cult_cols); N = len(cult_cols)
        fig, ax = plt.subplots(figsize=(7, 3.4))
        ax.hist(prev, bins=np.arange(0.5, N + 1.5, 1), color=FIG_BLUE, edgecolor="white")
        ax.set_xlabel("Number of cultivars in which a proteotype is present", fontsize=11)
        ax.set_ylabel("Number of proteotype clusters", fontsize=11)
        ax.set_xticks(range(1, N + 1))
        ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
        fig.savefig(out_png, bbox_inches="tight", dpi=600); plt.close(fig); return out_png
    except Exception as e:
        log.warning("prevalence hist failed: %s", e); return None


def plot_genomewide_per_cultivar(pav_df, cult_cols, out_png, log):
    try:
        M = pav_df[cult_cols].to_numpy(bool); prev = M.sum(1); N = len(cult_cols)
        core_mask = prev == N
        cc = [int((M[:, j] & core_mask).sum()) for j in range(N)]
        ac = [int((M[:, j] & ~core_mask).sum()) for j in range(N)]
        order = sorted(range(N), key=lambda j: -(cc[j] + ac[j])); x = np.arange(N)
        fig, ax = plt.subplots(figsize=(8, 3.6))
        ax.bar(x, [cc[j] for j in order], color=FIG_BLUE, label="core", edgecolor="white")
        ax.bar(x, [ac[j] for j in order], bottom=[cc[j] for j in order], color=FIG_AMBER, label="accessory", edgecolor="white")
        ax.set_xticks(x)
        ax.set_xticklabels([CULT_LABELS.get(cult_cols[j], cult_cols[j]) for j in order], rotation=45, ha="right", fontsize=8.5)
        ax.set_ylabel("Proteotype clusters", fontsize=10); ax.legend(frameon=False, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout(); fig.savefig(out_png, bbox_inches="tight", dpi=600); plt.close(fig); return out_png
    except Exception as e:
        log.warning("per-cultivar plot failed: %s", e); return None


def parse_uniprot_id(qseqid):
    parts = str(qseqid).split("|")
    if len(parts) >= 3:
        return parts[1], parts[2]
    return str(qseqid), ""


def build_proteotype_naming(inputs_dir, out_csv=None):
    pt = pd.read_csv(os.path.join(inputs_dir, "proteotypes.csv"))
    mem = pd.read_csv(os.path.join(inputs_dir, "members.csv"))
    hits = pd.read_csv(os.path.join(inputs_dir, "diamond_hits_all.csv"))
    hits = hits.sort_values(["evalue", "bitscore"], ascending=[True, False])
    best = hits.groupby(["cultivar", "sseqid"]).head(1).copy()
    parsed = [parse_uniprot_id(q) for q in best["qseqid"]]
    best["uacc"] = [a for a, _ in parsed]; best["uentry"] = [e for _, e in parsed]
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
    nm["composite_id"] = nm.apply(
        lambda r: f"{r.proteotype_id}__{r.representative_cultivar}__"
                  f"{r.uniprot_id if isinstance(r.uniprot_id, str) and r.uniprot_id else 'NA'}", axis=1)
    keep = ["proteotype_id", "family", "uniprot_id", "uniprot_entry", "representative_cultivar",
            "representative_source_id", "n_cultivars", "novel_vs_reference", "composite_id",
            "cultivars", "source_ids"]
    nm = nm[[c for c in keep if c in nm.columns]]
    if out_csv:
        nm.to_csv(out_csv, index=False)
    return nm


def _iter_fasta(path):
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header, chunks = None, []
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            elif header is not None:
                chunks.append(line.strip())
        if header is not None:
            yield header.split()[0], header, "".join(chunks)


def _write_fasta_record(fh, rec_id, seq, width=60):
    fh.write(f">{rec_id}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def _gene_id_from_header(header, rec_id):
    m = re.search(r"\bgene:(\S+)", header)
    if m:
        return m.group(1)
    return re.sub(r"\.\d+$", "", rec_id)


def _collapse_isoforms(in_fa, out_fa, log):
    best = {}
    for rid, header, seq in _iter_fasta(in_fa):
        g = _gene_id_from_header(header, rid)
        cur = best.get(g)
        if cur is None or len(seq) > cur[0] or (len(seq) == cur[0] and rid < cur[1]):
            best[g] = (len(seq), rid)
    keep = {rid for _, rid in best.values()}
    n_in = n_out = 0
    with open(out_fa, "w") as fo:
        for rid, _header, seq in _iter_fasta(in_fa):
            n_in += 1
            if rid in keep:
                keep.discard(rid)
                _write_fasta_record(fo, rid, seq)
                n_out += 1
    log.info("    isoform collapse: %d proteins -> %d genes (%s)",
             n_in, n_out, os.path.basename(out_fa))
    return out_fa, n_out


def _run_to_logfile(cmd, logfile, log, check=True):
    printable = " ".join(str(c) for c in cmd)
    log.debug("RUN: %s", printable)
    with open(logfile, "ab") as fh:
        fh.write(("\n$ " + printable + "\n").encode())
        fh.flush()
        proc = subprocess.run([str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT)
    if check and proc.returncode != 0:
        tail = ""
        try:
            with open(logfile, "rb") as fh2:
                fh2.seek(max(0, os.path.getsize(logfile) - 2000))
                tail = fh2.read().decode("utf-8", "replace")
        except Exception:
            pass
        hint = ""
        if proc.returncode == -9:
            hint = (" -- rc=-9 means the process was KILLED, almost always by the "
                    "out-of-memory killer. Lower whole_proteome.block_size_gb / "
                    "memory_limit_gb, or raise #SBATCH --mem.")
        log.error("Command failed (rc=%s)%s: %s\n%s", proc.returncode, hint, printable, tail)
        raise RuntimeError(f"Command failed: {cmd[0]} (rc={proc.returncode}){hint}")
    return proc


def _diamond_makedb(diamond, fasta, db, threads, log, logfile):
    _run_to_logfile([diamond, "makedb", "--in", fasta, "--db", db,
                     "--threads", str(threads), "--quiet"], logfile, log)
    return db + ".dmnd"


def _diamond_blastp_generic(diamond, query_fa, db, out_tsv, log, logfile, threads,
                            min_id, min_cov, evalue, max_target_seqs=1, sensitivity="",
                            block_size=None, index_chunks=None):
    cmd = [diamond, "blastp", "--query", query_fa, "--db", db + ".dmnd", "--out", out_tsv,
           "--outfmt", "6", *DIAMOND_FIELDS, "--evalue", str(evalue), "--id", str(min_id),
           "--query-cover", str(min_cov), "--max-target-seqs", str(max_target_seqs),
           "--max-hsps", "1", "--threads", str(threads), "--quiet"]
    if sensitivity:
        cmd += [s for s in str(sensitivity).split() if s]
    if block_size:
        cmd += ["--block-size", str(block_size)]
    if index_chunks:
        cmd += ["--index-chunks", str(index_chunks)]
    _run_to_logfile(cmd, logfile, log)
    return out_tsv


def _best_hits_from_tsv(tsv):
    best = {}
    if not (os.path.exists(tsv) and os.path.getsize(tsv) > 0):
        return {}
    with open(tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                bits = float(f[11])
            except ValueError:
                continue
            cur = best.get(f[0])
            if cur is None or bits > cur[0]:
                best[f[0]] = (bits, f[1])
    return {q: s for q, (_b, s) in best.items()}


class _UnionFind:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def _diamond_cluster(diamond, in_fa, out_tsv, log, logfile, threads,
                     approx_id, member_cover, memory_limit_gb):
    def _run(args, timeout=120):
        try:
            p = subprocess.run([str(diamond)] + args, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, timeout=timeout)
            return p, (p.stdout or "")
        except Exception as e:
            log.info("[whole-proteome] could not run 'diamond %s' (%s).",
                     " ".join(args), e)
            return None, ""

    _p, _cmds = _run(["help"])
    if _p is None:
        return None
    if "cluster" not in _cmds:
        log.warning("[whole-proteome] this DIAMOND build lists no 'cluster' command. "
                    "Falling back to blastp self-alignment + union-find, which is "
                    "TRANSITIVE and can chain dissimilar proteins into one cluster. "
                    "Treat the genome-wide partition as approximate.")
        return None
    probe, help_txt = _run(["cluster"])
    if probe is None:
        return None
    try:
        _v = subprocess.run([str(diamond), "--version"], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, timeout=60)
        log.info("[whole-proteome] %s", (_v.stdout or "").strip().splitlines()[0])
    except Exception:
        pass
    log.info("[whole-proteome] using 'diamond cluster' (greedy set-cover, NOT "
             "transitive) at approx-id %s / member-cover %s.", approx_id, member_cover)
    cmd = [diamond, "cluster", "-d", in_fa, "-o", out_tsv, "--threads", str(threads)]
    if "--approx-id" in help_txt:
        cmd += ["--approx-id", str(approx_id)]
    if "--member-cover" in help_txt:
        cmd += ["--member-cover", str(member_cover)]
    if memory_limit_gb and ("--memory-limit" in help_txt or "-M " in help_txt):
        cmd += ["--memory-limit", f"{memory_limit_gb}G"]
    cmd += ["--quiet"]
    if "--approx-id" not in help_txt and "--approx-id" not in cmd:
        cmd += ["--approx-id", str(approx_id), "--member-cover", str(member_cover)]
    proc = _run_to_logfile(cmd, logfile, log, check=False)
    if proc.returncode != 0:
        log.info("[whole-proteome] 'diamond cluster' returned rc=%s; "
                 "falling back to blastp clustering (see %s).",
                 proc.returncode, os.path.basename(logfile))
        return None
    if not (os.path.exists(out_tsv) and os.path.getsize(out_tsv) > 0):
        log.warning("[whole-proteome] 'diamond cluster' produced no output; "
                    "falling back to blastp clustering.")
        return None
    return out_tsv


def _cluster_novel_sequences(diamond, novel_fa, n_uniq, work, log, logfile, wp, threads):
    if n_uniq == 0:
        return []
    if not wp.get("cluster_novel", True):
        log.info("[whole-proteome] cluster_novel=false -> each unique novel protein "
                 "is its own cluster.")
        return [[i] for i in range(n_uniq)]

    uf = _UnionFind(n_uniq)
    linked = False

    ctsv = _diamond_cluster(diamond, novel_fa, os.path.join(work, "novel_clusters.tsv"),
                            log, logfile, threads,
                            wp.get("cluster_approx_id", 90), wp.get("cluster_member_cover", 80),
                            wp.get("memory_limit_gb", 48))
    if ctsv:
        n_edges = 0
        with open(ctsv) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                try:
                    a, b = int(f[0][1:]), int(f[1][1:])
                except ValueError:
                    continue
                if 0 <= a < n_uniq and 0 <= b < n_uniq:
                    uf.union(a, b)
                    n_edges += 1
        log.info("[whole-proteome] diamond cluster: %d membership edges.", n_edges)
        linked = n_edges > 0

    if not linked:
        try:
            log.warning("[whole-proteome] falling back to blastp self-alignment + "
                        "union-find for %d unique sequences. This is TRANSITIVE "
                        "clustering and is not equivalent to 'diamond cluster'; "
                        "chained merges inflate cluster sizes and the apparent core "
                        "fraction.", n_uniq)
            ndb = os.path.join(work, "novel_db")
            _diamond_makedb(diamond, novel_fa, ndb, threads, log, logfile)
            nout = os.path.join(work, "novel_self.tsv")
            _cid = float(wp.get("cluster_approx_id", 90))
            _ccov = float(wp.get("cluster_member_cover", 80))
            log.info("[whole-proteome] fallback clustering at %.0f%% identity / %.0f%% "
                     "coverage (cluster_approx_id, not min_identity).", _cid, _ccov)
            _diamond_blastp_generic(
                diamond, novel_fa, ndb, nout, log, logfile, threads,
                _cid, _ccov, wp.get("evalue", 1e-5),
                max_target_seqs=int(wp.get("fallback_max_target_seqs", 25)),
                sensitivity=wp.get("sensitivity", ""),
                block_size=wp.get("block_size_gb"), index_chunks=wp.get("index_chunks"))
            n_edges = 0
            with open(nout) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 2:
                        continue
                    try:
                        a, b = int(f[0][1:]), int(f[1][1:])
                    except ValueError:
                        continue
                    if a != b and 0 <= a < n_uniq and 0 <= b < n_uniq:
                        uf.union(a, b)
                        n_edges += 1
            log.info("[whole-proteome] blastp clustering: %d edges.", n_edges)
        except Exception as e:
            log.warning("[whole-proteome] novel clustering failed (%s); each unique novel "
                        "protein will be reported as its own cluster.", e)

    groups = {}
    for i in range(n_uniq):
        groups.setdefault(uf.find(i), []).append(i)
    out = list(groups.values())
    if out:
        sizes = sorted((len(g) for g in out), reverse=True)
        log.info("[whole-proteome] %d novel clusters; largest %d members, "
                 "top-5 sizes %s, singletons %d",
                 len(out), sizes[0], sizes[:5], sum(1 for s in sizes if s == 1))
        if sizes[0] > 0.05 * n_uniq:
            log.warning("[whole-proteome] the largest novel cluster holds %.1f%% of all "
                        "unique novel proteins. That is the signature of transitive "
                        "chaining through a permissive identity threshold. Check "
                        "whole_proteome.cluster_approx_id before using the genome-wide "
                        "partition.", 100.0 * sizes[0] / n_uniq)
    return out


def run_whole_proteome_reference(proteomes, ref_cult, out_dir, supp_dir, diamond, cfg, log):
    import hashlib

    wp = cfg.get("whole_proteome", {})
    if not wp.get("enabled", False):
        log.info("[whole-proteome] disabled in config; skipping genome-wide analysis.")
        return {}
    if ref_cult not in proteomes:
        log.warning("[whole-proteome] reference proteome '%s' unavailable; skipping.", ref_cult)
        return {}

    import json

    threads = cfg["homology"]["threads"]
    resume = bool(wp.get("resume", True))
    min_id = float(wp.get("min_identity", 40))
    min_cov = float(wp.get("min_coverage", 60))
    evalue = float(wp.get("evalue", 1e-5))
    work = wp.get("work_dir") or os.path.join(
        os.path.abspath(cfg["project"]["inputs_dir"]), "whole_proteome_work")
    os.makedirs(work, exist_ok=True)
    logfile = os.path.join(work, "diamond.log")

    cults = list(proteomes)
    cidx = {c: i for i, c in enumerate(cults)}

    sig = {"reference": ref_cult, "cultivars": sorted(cults),
           "collapse_isoforms": bool(wp.get("collapse_isoforms", True)),
           "min_identity": min_id, "min_coverage": min_cov, "evalue": evalue,
           "sensitivity": str(wp.get("sensitivity", ""))}
    sig_path = os.path.join(work, "settings.json")
    old_sig = None
    if os.path.exists(sig_path):
        try:
            with open(sig_path) as fh:
                old_sig = json.load(fh)
        except Exception:
            old_sig = None
    if resume and old_sig is not None and old_sig != sig:
        changed = sorted({k for k in set(old_sig) | set(sig) if old_sig.get(k) != sig.get(k)})
        log.warning("[whole-proteome] cached work in %s was produced with different "
                    "settings (%s changed); ignoring the cache and recomputing.",
                    work, ", ".join(changed))
        resume = False
    with open(sig_path, "w") as fh:
        json.dump(sig, fh, indent=2, sort_keys=True)
    log.info("[whole-proteome] work dir: %s (resume=%s)", work, resume)

    seqs = dict(proteomes)
    if wp.get("collapse_isoforms", True):
        for c in cults:
            tgt = os.path.join(work, f"{sanitize_name(c)}.longest.faa")
            if resume and os.path.exists(tgt) and os.path.getsize(tgt) > 0:
                log.info("    reuse collapsed proteome: %s", os.path.basename(tgt))
            else:
                _collapse_isoforms(proteomes[c], tgt, log)
            seqs[c] = tgt

    refdb = os.path.join(work, "ref_db")
    if not (resume and os.path.exists(refdb + ".dmnd")):
        _diamond_makedb(diamond, seqs[ref_cult], refdb, threads, log, logfile)

    ref_ids = [rid for rid, _h, _s in _iter_fasta(seqs[ref_cult])]
    ridx = {rid: i for i, rid in enumerate(ref_ids)}
    log.info("[whole-proteome] reference '%s': %d proteins.", ref_cult, len(ref_ids))

    M = np.zeros((len(ref_ids), len(cults)), dtype=bool)
    M[:, cidx[ref_cult]] = True

    novel_fa = os.path.join(work, "novel_unique.faa")
    seen = {}
    novel_mask = []
    n_novel_total = 0
    with open(novel_fa, "w") as fo:
        for cult in cults:
            if cult == ref_cult:
                continue
            tsv = os.path.join(work, f"{sanitize_name(cult)}_vs_ref.tsv")
            if resume and os.path.exists(tsv) and os.path.getsize(tsv) > 0:
                log.info("    reuse vs-reference hits: %s", os.path.basename(tsv))
            else:
                _diamond_blastp_generic(
                    diamond, seqs[cult], refdb, tsv, log, logfile, threads,
                    min_id, min_cov, evalue, max_target_seqs=1,
                    sensitivity=wp.get("sensitivity", ""),
                    block_size=wp.get("block_size_gb"), index_chunks=wp.get("index_chunks"))
            best = _best_hits_from_tsv(tsv)
            n_novel = 0
            for rid, _h, seq in _iter_fasta(seqs[cult]):
                tgt = best.get(rid)
                if tgt is not None:
                    j = ridx.get(tgt)
                    if j is not None:
                        M[j, cidx[cult]] = True
                    continue
                if not seq:
                    continue
                h = hashlib.md5(seq.encode()).digest()
                u = seen.get(h)
                if u is None:
                    u = len(novel_mask)
                    seen[h] = u
                    novel_mask.append(0)
                    _write_fasta_record(fo, f"U{u}", seq)
                novel_mask[u] |= (1 << cidx[cult])
                n_novel += 1
            n_novel_total += n_novel
            log.info("  %-14s mapped=%d novel=%d", cult, len(best), n_novel)
            del best
    log.info("[whole-proteome] novel proteins: %d total -> %d unique sequences "
             "(exact duplicates merged).", n_novel_total, len(novel_mask))

    groups = _cluster_novel_sequences(diamond, novel_fa, len(novel_mask), work,
                                      log, logfile, wp, threads)

    if groups:
        NM = np.zeros((len(groups), len(cults)), dtype=bool)
        for k, members in enumerate(groups):
            mask = 0
            for i in members:
                mask |= novel_mask[i]
            for j in range(len(cults)):
                if mask & (1 << j):
                    NM[k, j] = True
        full = np.vstack([M, NM])
    else:
        full = M
    ids = ref_ids + [f"NOVEL_{k + 1:06d}" for k in range(len(groups))]
    gw = pd.DataFrame(full.astype(np.int8), columns=cults)
    gw.insert(0, "cluster_id", ids)
    gw["novel_vs_reference"] = [False] * len(ref_ids) + [True] * len(groups)
    gw_path = os.path.join(supp_dir, "genome_wide_pav.csv")
    gw.to_csv(gw_path, index=False)

    part, _prev = _pangenome_partition(gw, cults)
    part_path = os.path.join(supp_dir, "pangenome_partition.csv")
    pd.DataFrame({"category": list(part), "n_clusters": list(part.values())}).to_csv(
        part_path, index=False)
    log.info("[whole-proteome] %d reference + %d novel clusters; partition=%s",
             len(ref_ids), len(groups), part)

    arts = [("Genome-wide PAV (whole CS reference)", gw_path),
            ("Pan-genome partition table", part_path)]

    if wp.get("write_catalog", True):
        catalog = os.path.join(out_dir, "WholeProteome_catalog.annotated.faa")
        N = len(cults)
        soft = math.ceil(0.9 * N)

        def _cls(p):
            return ("core" if p == N else
                    "soft_core" if p >= soft else
                    "shell" if p >= 2 else "cloud")

        n_written = 0
        try:
            with open(catalog, "w") as fo:
                for rid, _hdr, seq in _iter_fasta(seqs[ref_cult]):
                    row = ridx.get(rid)
                    if row is None or not seq:
                        continue
                    prev = int(full[row].sum())
                    hdr = (f"{rid}__{ref_cult}__CSref family=whole_proteome "
                           f"rep_cultivar={ref_cult} n_cultivars={prev} "
                           f"prevalence_class={_cls(prev)} novel=False reference_locus={rid}")
                    _write_fasta_record(fo, hdr, seq)
                    n_written += 1
                if groups:
                    rep_of = {}
                    for k, members in enumerate(groups):
                        if members:
                            rep_of[members[0]] = k
                    for rid, _hdr, seq in _iter_fasta(novel_fa):
                        if not rid.startswith("U") or not seq:
                            continue
                        try:
                            u = int(rid[1:])
                        except ValueError:
                            continue
                        k = rep_of.get(u)
                        if k is None:
                            continue
                        prev = int(NM[k].sum())
                        rep_c = next((cults[j] for j in range(N) if NM[k, j]), cults[0])
                        cid = f"NOVEL_{k + 1:06d}"
                        hdr = (f"{cid}__{rep_c}__novel family=whole_proteome "
                               f"rep_cultivar={rep_c} n_cultivars={prev} "
                               f"prevalence_class={_cls(prev)} novel=True reference_locus=NA")
                        _write_fasta_record(fo, hdr, seq)
                        n_written += 1
            log.info("[whole-proteome] wrote annotated catalogue: %d cluster records -> %s",
                     n_written, os.path.basename(catalog))
            arts.append(("Whole-proteome annotated FASTA (all clusters)", catalog))
        except Exception as e:
            log.warning("[whole-proteome] catalogue FASTA failed: %s", e)

    pth = plot_pangenome_overview(gw, cults, os.path.join(out_dir, "Figure_pangenome_overview.png"), log)
    if pth:
        arts.append(("Figure: whole-proteome pan-genome overview", pth))
    pth = plot_prevalence_hist(gw, cults, os.path.join(supp_dir, "genomewide_prevalence_spectrum.png"), log)
    if pth:
        arts.append(("Suppl. Figure: genome-wide prevalence spectrum", pth))
    pth = plot_genomewide_per_cultivar(gw, cults, os.path.join(supp_dir, "genomewide_per_cultivar.png"), log)
    if pth:
        arts.append(("Suppl. Figure: genome-wide per-cultivar composition", pth))
    return {"genome_wide_pav": gw_path, "partition": part, "artifacts": arts}


def run_pipeline(cfg, args, log):
    ts = get_timestamp()
    outputs_root = os.path.abspath(cfg["project"]["outputs_dir"])
    inputs_dir = os.path.abspath(cfg["project"]["inputs_dir"])
    out_dir = os.path.join(outputs_root, f"pan_proteome_build_{ts}")
    os.makedirs(out_dir, exist_ok=True)
    supp_dir = os.path.join(out_dir, "supplementary")
    os.makedirs(supp_dir, exist_ok=True)
    os.makedirs(inputs_dir, exist_ok=True)

    if not _HAVE_BIO:
        raise RuntimeError("Biopython is required for the full pipeline. "
                           "Activate the conda environment first.")
    _install_proxy(cfg, log)

    diamond = require_tool("diamond", log)
    mafft = require_tool("mafft", log)
    fasttree = shutil.which("FastTree") or shutil.which("fasttree")
    if not fasttree:
        log.warning("FastTree not found; trees will be skipped.")

    log.info("Tool versions: diamond=%s | mafft=%s | fasttree=%s",
             tool_version("diamond", ["--version"], log),
             tool_version("mafft", ["--version"], log),
             tool_version(os.path.basename(fasttree) if fasttree else "fasttree", [], log))

    only = None
    if args.cultivars:
        only = {c.strip().lower() for c in args.cultivars.split(",")}
        cfg["cultivars"] = {k: v for k, v in cfg["cultivars"].items() if k in only}

    log.info("=== [1/7] Reference sequences ===")
    refs = fetch_reference_sequences(cfg, os.path.join(inputs_dir, "references"), log)

    log.info("=== [2/7] Cultivar proteomes ===")
    base_url, release_dir = resolve_release(cfg, log)
    proteomes = fetch_cultivar_proteomes(cfg, inputs_dir, base_url, release_dir, log, only=only)
    cultivars = list(proteomes)

    log.info("=== [3/7] Homology search + extraction ===")
    all_members = []
    hits_frames = []
    dwork = os.path.join(out_dir, "diamond")
    for cult in cultivars:
        wd = os.path.join(dwork, sanitize_name(cult))
        hits = diamond_search(diamond, refs["combined"], proteomes[cult], wd, cfg, log)
        if not hits.empty:
            h = hits.copy()
            h.insert(0, "cultivar", cult)
            hits_frames.append(h)
        fam_map = assign_families(hits, refs["qid2family"])
        members = extract_members(cult, proteomes[cult], fam_map, log)
        log.info("  %-14s hits=%d  RCA=%d  CPN60=%d", cult, len(fam_map),
                 sum(1 for m in members if m["family"] == "RCA"),
                 sum(1 for m in members if m["family"] == "CPN60"))
        all_members.extend(members)
    if not all_members:
        raise RuntimeError("No RCA/Cpn60 homologs found in any proteome. "
                           "Relax homology thresholds or check reference queries.")

    log.info("=== [4/7] Non-redundant proteotypes + PAV ===")
    ref_cult = "chinese_spring" if "chinese_spring" in cultivars else cultivars[0]
    if ref_cult != "chinese_spring":
        log.warning("Chinese Spring not in panel; using '%s' as novelty baseline.", ref_cult)
    proteotypes, member_records = dedupe_to_proteotypes(all_members, reference_cultivar=ref_cult)
    pav_df = build_pav_matrix(proteotypes, cultivars)

    ext_cache = (cfg.get("external", {}).get("cache_dir")
                 or os.path.join(inputs_dir, "external_cache"))
    os.makedirs(ext_cache, exist_ok=True)

    log.info("=== [5/7] Alignment + tree per family ===")
    aln_by_family = {}
    tree_support = {}
    artifacts = []
    rooting_report = []
    for fam in FAMILIES:
        fam_fa = os.path.join(out_dir, f"{fam}_proteotypes.faa")
        n = write_family_fasta(proteotypes, fam, fam_fa)
        if n and cfg.get("external", {}).get("fetch_outgroups", False):
            try:
                import external_data
                og_ids = external_data.append_outgroups(
                    fam_fa, cfg["external"]["outgroups"].get(fam, []),
                    ext_cache, cfg, log)
                if og_ids:
                    n += len(og_ids)
                    cfg["alignment"].setdefault("outgroup", {})[fam] = og_ids
                    log.info("  %s: %d outgroup sequence(s) appended (%s)",
                             fam, len(og_ids), ", ".join(og_ids))
                    cfg["alignment"]["root"] = "outgroup"
            except Exception as e:
                log.warning("outgroup fetch failed for %s (%s); tree will be "
                            "midpoint-rooted instead.", fam, e)
        if n == 0:
            continue
        artifacts.append((f"{fam} proteotype sequences", fam_fa))
        if n >= 2:
            aln = os.path.join(out_dir, f"{fam}.aln.fasta")
            run_mafft(mafft, fam_fa, aln, cfg, log)
            artifacts.append((f"{fam} alignment", aln))
            recs = {r.id: str(r.seq) for r in SeqIO.parse(aln, "fasta")}
            aln_by_family[fam] = recs
            if fasttree and n >= 3:
                nwk = os.path.join(out_dir, f"{fam}.nwk")
                run_fasttree(fasttree, aln, nwk, log)
                artifacts.append((f"{fam} tree (Newick)", nwk))
                novel_ids = {p["proteotype_id"] for p in proteotypes
                             if p.get("family") == fam and p.get("novel_vs_reference")}
                png = os.path.join(out_dir, f"{fam}_tree_circular.png")
                try:
                    import tree_tools
                    troot = tree_tools.parse_newick(open(nwk).read())
                    summ = tree_tools.support_summary(troot)
                    if summ.get("n"):
                        log.info("  %s tree support: median %.3f; %d/%d nodes >= 0.80",
                                 fam, summ["median"], summ["ge_0.80"], summ["n"])
                        tree_support[fam] = summ
                    thr = float(cfg["alignment"].get("collapse_support_below", 0.8))
                    n_col = tree_tools.collapse_unsupported(troot, thr)
                    mode = cfg["alignment"].get("root", "midpoint")
                    tipnames = {l["name"] for l in tree_tools.leaves(troot)}
                    og = [t for t in (cfg["alignment"].get("outgroup", {}).get(fam) or [])
                          if t in tipnames]
                    want = mode
                    if mode == "outgroup" and og:
                        troot = tree_tools.outgroup_root(troot, og)
                    elif mode in ("outgroup", "midpoint"):
                        troot = tree_tools.midpoint_root(troot)
                        mode = "midpoint"
                    if want == "outgroup" and mode != "outgroup":
                        cfgd = cfg["alignment"].get("outgroup", {}).get(fam) or []
                        why = ("no outgroup accessions configured" if not cfgd else
                               "configured outgroup ids %s are not tip names in the tree "
                               "(tips look like %s)" % (
                                   cfgd, sorted(tipnames)[:3]))
                        log.error("  *** %s: outgroup rooting REQUESTED but NOT achieved "
                                  "-> tree is MIDPOINT-rooted. Reason: %s", fam, why)
                        log.error("  *** Outgroup rooting was not achieved for %s.", fam)
                    rooting_report.append({"family": fam, "requested": want,
                                           "achieved": mode,
                                           "outgroup_tips": ";".join(og),
                                           "n_collapsed": n_col})
                    sup_nwk = os.path.join(out_dir, f"{fam}.supported.nwk")
                    with open(sup_nwk, "w") as fh:
                        fh.write(tree_tools.to_newick(troot) + ";\n")
                    artifacts.append((f"{fam} tree, collapsed and rooted (Newick)", sup_nwk))
                    tree_tools.render_circular(
                        troot, png, novel_ids=novel_ids, prefix=(fam + "_"),
                        title_note=(f"Node labels: SH-like local support; nodes < {thr:g} "
                                    f"collapsed ({n_col}); {mode}-rooted"))
                except Exception as e:
                    log.warning("tree_tools rendering failed (%s); using the plain renderer.", e)
                    png = plot_circular_tree(nwk, png,
                                             f"{FAM_DISP.get(fam, fam)} proteotype tree", log,
                                             novel_ids=novel_ids, family=fam)
                if png:
                    artifacts.append((f"Figure: {fam} tree", png))

    if rooting_report:
        _rr = os.path.join(out_dir, "tree_rooting.csv")
        pd.DataFrame(rooting_report).to_csv(_rr, index=False)
        artifacts.append(("Tree rooting report", _rr))
        for _r in rooting_report:
            log.info("  rooting: %s requested=%s achieved=%s outgroups=[%s]",
                     _r["family"], _r["requested"], _r["achieved"],
                     _r["outgroup_tips"])
        _bad = [r["family"] for r in rooting_report
                if r["requested"] == "outgroup" and r["achieved"] != "outgroup"]
        if _bad:
            log.error("  *** %d famil(y/ies) fell back to midpoint rooting: %s",
                      len(_bad), ", ".join(_bad))
            log.error("  *** See tree_rooting.csv. Fix before quoting outgroup rooting.")

    log.info("=== [6/7] Descriptors + statistics ===")
    feat_rows = []
    for pt in proteotypes:
        f = physicochemical_features(pt["seq"])
        f.update({"proteotype_id": pt["proteotype_id"], "family": pt["family"],
                  "novel_vs_reference": pt["novel_vs_reference"]})
        feat_rows.append(f)
    features_df = pd.DataFrame(feat_rows)
    pt2feat = {r["proteotype_id"]: r for r in feat_rows}
    member_feat_rows = []
    for m in member_records:
        base = dict(pt2feat.get(m["proteotype_id"], {}))
        base.update({"cultivar": m["cultivar"], "family": m["family"],
                     "source_id": m["source_id"], "proteotype_id": m["proteotype_id"]})
        member_feat_rows.append(base)
    member_features_df = pd.DataFrame(member_feat_rows)

    resistant, sensitive = groups_from_config(cfg)
    resistant = [c for c in resistant if c in cultivars]
    sensitive = [c for c in sensitive if c in cultivars]
    group_stats = compare_group_features(member_features_df, resistant, sensitive, log) \
        if (resistant and sensitive) else pd.DataFrame(columns=GROUP_COMPARISON_COLS)

    poly_frames = []
    pos_assoc = {}
    for fam, recs in aln_by_family.items():
        poly_frames.append(polymorphic_sites(recs, fam))
        if resistant and sensitive:
            pos_assoc[fam] = position_group_association(recs, member_records,
                                                        resistant, sensitive, fam, log)
    poly_df = pd.concat(poly_frames, ignore_index=True) if poly_frames \
        else pd.DataFrame(columns=POLY_COLS)
    known_df = map_known_residues(cfg, aln_by_family, member_records, log)

    log.info("=== [7/7] Persistence, plots, report ===")
    catalog = os.path.join(out_dir, "PanProteome_catalog.faa")
    SeqIO.write([SeqRecord(Seq(p["seq"]), id=p["proteotype_id"],
                           description=f"{p['family']} novel={p['novel_vs_reference']} "
                                       f"cultivars={p['cultivars']}") for p in proteotypes],
                catalog, "fasta")
    artifacts.insert(0, ("Pan-proteome catalog (all proteotypes)", catalog))

    proteotypes_df = pd.DataFrame([{k: v for k, v in p.items() if k != "seq"}
                                   for p in proteotypes])
    members_df = pd.DataFrame([{k: v for k, v in m.items() if k != "seq"}
                               for m in member_records])
    csv_specs = [
        ("proteotypes.csv", proteotypes_df), ("members.csv", members_df),
        ("pav_matrix.csv", pav_df), ("physicochemical_features.csv", features_df),
        ("polymorphic_sites.csv", poly_df), ("group_comparison.csv", group_stats),
        ("known_residue_overlay.csv", known_df),
    ]
    for fam, pa in pos_assoc.items():
        csv_specs.append((f"position_association_{fam}.csv", pa))
    for fname, df in csv_specs:
        if df is not None:
            path = os.path.join(out_dir, fname)
            df.to_csv(path, index=False)
            artifacts.append((fname.replace(".csv", "").replace("_", " ").title(), path))

    ov = plot_overview(proteotypes, pav_df, cultivars,
                       os.path.join(out_dir, "Figure_overview.png"), log)
    if ov:
        artifacts.append(("Figure: pan-proteome overview", ov))
    for fam in FAMILIES:
        pn = plot_novel_pav(proteotypes, pav_df, fam, cultivars,
                            os.path.join(out_dir, f"{fam}_novel_PAV.png"), log)
        if pn:
            artifacts.append((f"Figure: {fam} novel PAV", pn))
        p1 = plot_pav_heatmap(pav_df, fam, os.path.join(supp_dir, f"{fam}_PAV_heatmap.png"), log)
        if p1:
            artifacts.append((f"Suppl. Figure: {fam} PAV heatmap", p1))
        if cfg.get("report", {}).get("draw_circos", False):
            p2 = plot_circos_pav(pav_df, fam, os.path.join(supp_dir, f"{fam}_PAV_circos.png"), log)
            if p2:
                artifacts.append((f"Suppl. Figure: {fam} circos PAV", p2))
        if not poly_df.empty:
            p3 = plot_entropy(poly_df, fam, os.path.join(out_dir, f"{fam}_variability.png"), log)
            if p3:
                artifacts.append((f"Figure: {fam} variability", p3))
        if resistant and sensitive:
            p4 = plot_feature_boxplots(member_features_df, resistant, sensitive, fam,
                                       os.path.join(supp_dir, f"{fam}_descriptors_by_group.png"), log)
            if p4:
                artifacts.append((f"Figure: {fam} descriptors by group", p4))

    db_path = os.path.join(out_dir, "wheat_pan_proteome.db")
    run_meta = pd.DataFrame([{"timestamp": ts, "release": release_dir,
                              "n_cultivars": len(cultivars),
                              "cultivars": ";".join(cultivars),
                              "resistant": ";".join(resistant),
                              "sensitive": ";".join(sensitive),
                              "n_proteotypes": len(proteotypes)}])
    write_sqlite(db_path, {"Proteotypes": proteotypes_df, "Members": members_df,
                           "PAV": pav_df, "Features": features_df,
                           "Polymorphic": poly_df, "GroupComparison": group_stats,
                           "KnownResidues": known_df, "RunMeta": run_meta}, log)
    artifacts.append(("SQLite catalog", db_path))
    if hits_frames:
        allhits = os.path.join(out_dir, "diamond_hits_all.csv")
        pd.concat(hits_frames, ignore_index=True).to_csv(allhits, index=False)
        artifacts.append(("All DIAMOND hits", allhits))

    try:
        naming_csv = os.path.join(out_dir, "proteotype_naming.csv")
        naming_df = build_proteotype_naming(out_dir, naming_csv)
        shutil.copy(naming_csv, os.path.join(supp_dir, "Table_S7_proteotype_naming.csv"))
        artifacts.append(("Proteotype naming table (UniProt+cultivar+id)", naming_csv))
        info = {r.proteotype_id: r for r in naming_df.itertuples()}
        ann = os.path.join(out_dir, "PanProteome_catalog.annotated.faa")
        with open(ann, "w") as fo:
            for pp in proteotypes:
                r = info.get(pp["proteotype_id"])
                if r is None:
                    cid, desc = pp["proteotype_id"], f"family={pp['family']}"
                else:
                    cid = r.composite_id
                    desc = (f"family={r.family} uniprot={r.uniprot_id} uniprot_entry={r.uniprot_entry} "
                            f"rep_cultivar={r.representative_cultivar} rep_source_id={r.representative_source_id} "
                            f"n_cultivars={r.n_cultivars} novel={r.novel_vs_reference}")
                seq = pp["seq"]
                fo.write(f">{cid} {desc}\n")
                for i in range(0, len(seq), 60):
                    fo.write(seq[i:i + 60] + "\n")
        artifacts.append(("Annotated proteotype FASTA (composite ids)", ann))
    except Exception as e:
        log.warning("proteotype naming / annotated FASTA failed: %s", e)

    for _f in ("members.csv", "diamond_hits_all.csv", "physicochemical_features.csv",
               "polymorphic_sites.csv"):
        _src = os.path.join(out_dir, _f)
        if os.path.exists(_src):
            shutil.copy(_src, os.path.join(supp_dir, _f))

    try:
        gw_res = run_whole_proteome_reference(proteomes, ref_cult, out_dir, supp_dir, diamond, cfg, log)
        for _label, _path in gw_res.get("artifacts", []):
            artifacts.append((_label, _path))
    except Exception as e:
        log.warning("whole-proteome analysis failed: %s", e)

    audit_out = {}
    try:
        import panproteome_audit
        audit_dir = os.path.join(out_dir, "audit")
        audit_out = panproteome_audit.run_audit(
            out_dir, audit_dir, reference_cultivar=ref_cult,
            group2_clade=(cfg.get("report", {}).get("candidate_clade") or None))
        for fn, label in [
                ("locus_counts_by_cultivar.csv", "Protein models and gene loci per cultivar"),
                ("cs_locus_table.csv", "Reference gene-locus table"),
                ("specificity_audit.csv", "Recruitment specificity audit"),
                ("cpn60_classes.csv", "Cpn60 provisional class assignment"),
                ("novel_divergence.csv", "Divergence of reference-absent proteotypes"),
                ("granularity_partition.csv", "Prevalence partition at three granularities"),
                ("core_fraction_ci.csv", "Bootstrap CI on core fractions"),
                ("candidate_table.csv", "Prioritised candidate proteotypes")]:
            p = os.path.join(audit_dir, fn)
            if os.path.exists(p):
                artifacts.append((label, p))
    except Exception as e:
        log.warning("post-hoc audit failed: %s", e)

    try:
        if cfg.get("external", {}).get("annotate_proteotypes", False):
            import external_data
            ann = external_data.annotate_naming_table(
                os.path.join(out_dir, "proteotype_naming.csv"), ext_cache, cfg, log)
            if ann is not None:
                p = os.path.join(out_dir, "proteotype_annotation.csv")
                ann.to_csv(p, index=False)
                shutil.copy(p, os.path.join(supp_dir, "proteotype_annotation.csv"))
                artifacts.append(("UniProt/InterPro annotation of proteotypes", p))
                audit_out["annotation"] = ann
    except Exception as e:
        log.warning("external annotation failed: %s", e)

    print("\n" + "=" * 70)
    print(" PAN-PROTEOME BUILD COMPLETE")
    print("=" * 70)
    print(f" Output directory : {out_dir}")
    for fam in FAMILIES:
        fp = [p for p in proteotypes if p["family"] == fam]
        print(f" {fam:6s}: {len(fp)} proteotypes "
              f"({sum(1 for p in fp if p['novel_vs_reference'])} novel vs baseline)")
    if not (resistant and sensitive):
        print(" NOTE: set resistant/sensitive groups in config.yaml to enable the comparison.")
    print("=" * 70)
    return out_dir


def selftest(log):
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    poly_a = "A" * 100
    check("aliphatic_index(poly-Ala)==100", abs(aliphatic_index(poly_a) - 100.0) < 1e-6)
    check("gravy(poly-Ala)==1.8", abs(gravy(poly_a) - 1.8) < 1e-6)
    check("aromaticity(FYW mix)>0", aromaticity("FYWAAAA") > 0)
    check("mw positive", molecular_weight("ACDEFGHIKLMNPQRSTVWY") > 2000)
    check("pI acidic<7", isoelectric_point("DDDEEE") < 7)
    check("pI basic>7", isoelectric_point("KKKRRR") > 7)
    feats = physicochemical_features("MSTAVLGKRDDEEPPGGWWYY*X")
    check("features has aliphatic_index", "aliphatic_index" in feats)
    check("features length_std_aa < length", feats["length_std_aa"] <= feats["length"])

    members = [
        {"cultivar": "chinese_spring", "family": "RCA", "source_id": "cs1",
         "pident_to_ref": 100, "evalue_to_ref": 0.0, "length": 5, "seq": "MSTAA"},
        {"cultivar": "jagger", "family": "RCA", "source_id": "j1",
         "pident_to_ref": 99, "evalue_to_ref": 0.0, "length": 5, "seq": "MSTAA"},
        {"cultivar": "mace", "family": "RCA", "source_id": "m1",
         "pident_to_ref": 98, "evalue_to_ref": 0.0, "length": 5, "seq": "MSTAV"},
        {"cultivar": "mace", "family": "CPN60", "source_id": "m2",
         "pident_to_ref": 97, "evalue_to_ref": 0.0, "length": 6, "seq": "GGGGGG"},
    ]
    pts, recs = dedupe_to_proteotypes(members, reference_cultivar="chinese_spring")
    rca = [p for p in pts if p["family"] == "RCA"]
    check("RCA collapses to 2 proteotypes", len(rca) == 2)
    novel = [p for p in rca if p["novel_vs_reference"]]
    check("MSTAV flagged novel vs CS", any(p["cultivars"] == "mace" for p in novel))
    pav = build_pav_matrix(pts, ["chinese_spring", "jagger", "mace"])
    check("PAV matrix shape", pav.shape[0] == len(pts) and "mace" in pav.columns)

    aln = {"a": "MSTAA", "b": "MSTAA", "c": "MSTAV"}
    poly = polymorphic_sites(aln, "RCA")
    check("one polymorphic column at pos 5", len(poly) == 1 and int(poly.iloc[0]["aln_column"]) == 5)

    q = benjamini_hochberg([0.01, 0.02, 0.03, 0.9])
    check("BH monotone & <=1", (q <= 1).all() and q[-1] >= q[0])

    if _HAVE_SCIPY:
        mf = pd.DataFrame([
            {"cultivar": "r1", "family": "RCA", "aliphatic_index": 90, "gravy": 0.1,
             "instability_index": 30, "frac_charged_DEKRH": 0.2, "net_charge_frac": 0.01,
             "pI": 6.0, "frac_acidic_DE": 0.1, "frac_basic_KR": 0.11, "frac_proline": 0.05,
             "frac_glycine": 0.06, "frac_arg_of_RK": 0.5, "aromaticity": 0.08},
            {"cultivar": "r2", "family": "RCA", "aliphatic_index": 92, "gravy": 0.12,
             "instability_index": 31, "frac_charged_DEKRH": 0.21, "net_charge_frac": 0.02,
             "pI": 6.1, "frac_acidic_DE": 0.1, "frac_basic_KR": 0.12, "frac_proline": 0.05,
             "frac_glycine": 0.06, "frac_arg_of_RK": 0.5, "aromaticity": 0.08},
            {"cultivar": "s1", "family": "RCA", "aliphatic_index": 70, "gravy": -0.2,
             "instability_index": 45, "frac_charged_DEKRH": 0.3, "net_charge_frac": -0.02,
             "pI": 5.4, "frac_acidic_DE": 0.16, "frac_basic_KR": 0.14, "frac_proline": 0.05,
             "frac_glycine": 0.06, "frac_arg_of_RK": 0.5, "aromaticity": 0.08},
            {"cultivar": "s2", "family": "RCA", "aliphatic_index": 71, "gravy": -0.18,
             "instability_index": 44, "frac_charged_DEKRH": 0.31, "net_charge_frac": -0.03,
             "pI": 5.3, "frac_acidic_DE": 0.17, "frac_basic_KR": 0.14, "frac_proline": 0.05,
             "frac_glycine": 0.06, "frac_arg_of_RK": 0.5, "aromaticity": 0.08},
        ])
        gc = compare_group_features(mf, ["r1", "r2"], ["s1", "s2"], log)
        check("group comparison returns rows", not gc.empty)
        check("group comparison has q_value_BH", "q_value_BH" in gc.columns)

    print("\nSELFTEST:", "ALL PASSED" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


def check_env(cfg, log):
    print("Environment check")
    print("  python           :", sys.version.split()[0])
    for mod, flag in [("Biopython", _HAVE_BIO), ("scipy", _HAVE_SCIPY),
                      ("seaborn", _HAVE_SNS), ("ProtParam(instability)", _HAVE_PROTPARAM)]:
        print(f"  {mod:22s}: {'yes' if flag else 'NO'}")
    for tool in ["diamond", "mafft", "FastTree", "fasttree"]:
        print(f"  tool {tool:16s}: {shutil.which(tool) or 'not found'}")
    r, s = groups_from_config(cfg)
    print(f"  cultivars in config : {len(cfg['cultivars'])}")
    print(f"  resistant / sensitive: {len(r)} / {len(s)}",
          "(comparison ENABLED)" if (len(r) >= 2 and len(s) >= 2) else "(comparison disabled)")
    return 0


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(
        description="Extended pan-proteome of RCA and Cpn60 across bread wheat cultivars.")
    ap.add_argument("--config", default=os.path.join(here, "config.yaml"),
                    help="Path to config.yaml (default: alongside this script).")
    ap.add_argument("--cultivars", help="Comma-separated subset of cultivar keys to run.")
    ap.add_argument("--inputs", help="Override project.inputs_dir (persistent cache).")
    ap.add_argument("--outputs", help="Override project.outputs_dir.")
    ap.add_argument("--threads", type=int, help="Override homology.threads.")
    ap.add_argument("--selftest", action="store_true",
                    help="Run offline logic tests and exit (no network/tools needed).")
    ap.add_argument("--check-env", action="store_true",
                    help="Report tool/library availability and config status, then exit.")
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest(setup_logging(os.path.join(here, "outputs")))

    cfg = load_config(args.config)
    if args.inputs:
        cfg["project"]["inputs_dir"] = args.inputs
    if args.outputs:
        cfg["project"]["outputs_dir"] = args.outputs
    if args.threads:
        cfg["homology"]["threads"] = args.threads

    log = setup_logging(os.path.abspath(cfg["project"]["outputs_dir"]))
    if args.check_env:
        return check_env(cfg, log)

    try:
        run_pipeline(cfg, args, log)
        return 0
    except Exception as e:
        log.exception("Pipeline failed: %s", e)
        return 1


if __name__ == "__main__":
    sys.exit(main())
