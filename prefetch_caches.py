#!/usr/bin/env python3
"""Download everything the compute node cannot fetch itself. LOGIN NODE ONLY.

Setonix compute nodes have no internet. domain_profiling.py needs Pfam HMMs and
Swiss-Prot seed sequences; transcript_crossref.py needs the White et al. (2025) Zenodo
tables. None of those depend on the pipeline's own output, so they can be fetched before
the job is submitted and used offline inside it.

Usage
    python prefetch_caches.py                       # caches into ./inputs
    python prefetch_caches.py --inputs inputs --expression
"""
import argparse, os, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="inputs")
    ap.add_argument("--expression", action="store_true",
                    help="also fetch the 273 MB DESeq2 counts table")
    ap.add_argument("--record", default=None)
    a = ap.parse_args()
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    hmm = os.path.abspath(os.path.join(a.inputs, "hmm_cache"))
    zen = os.path.abspath(os.path.join(a.inputs, "zenodo_cache"))
    os.makedirs(hmm, exist_ok=True); os.makedirs(zen, exist_ok=True)
    rc = 0

    print("=" * 70); print(" 1/2  Pfam HMMs + Swiss-Prot class seeds ->", hmm); print("=" * 70)
    try:
        import domain_profiling as D
        D.get_pfam_hmms(hmm)
        seeds = D.get_seed_sets(hmm)
        for k, v in sorted(seeds.items()):
            n = v.count(">") if isinstance(v, str) else len(v)
            print(f"    {k:<24} {n} seed sequence(s)")
        missing = [k for k, v in seeds.items()
                   if (v.count('>') if isinstance(v, str) else len(v)) < D.MIN_SEEDS]
        if missing:
            print(f"    WARNING under {D.MIN_SEEDS} seeds for: {missing} "
                  "-> those class profiles will be skipped")
    except Exception as e:
        print(f"    FAILED: {type(e).__name__}: {e}"); rc = 1

    print(); print("=" * 70); print(" 2/2  White et al. (2025) Zenodo tables ->", zen); print("=" * 70)
    try:
        import transcript_crossref as T
        rec = a.record or T.ZENODO_RECORD
        avail = T.zenodo_files(rec)
        want = list(T.CORE_FILES) + (T.EXPR_FILES if a.expression else [])
        for name in want:
            if name not in avail:
                print(f"    WARNING {name} not in record {rec}"); continue
            url, size = avail[name]
            print(f"    {name} ({size/1e6:.1f} MB)")
            T.fetch(url, os.path.join(zen, name))
    except Exception as e:
        print(f"    FAILED: {type(e).__name__}: {e}"); rc = 1

    print(); print("=" * 70)
    for d, label in ((hmm, "hmm_cache"), (zen, "zenodo_cache")):
        fs = [f for f in os.listdir(d) if not f.endswith(".part")]
        tot = sum(os.path.getsize(os.path.join(d, f)) for f in fs)
        print(f"  {label:<14} {len(fs):>3} file(s)  {tot/1e6:>8.1f} MB")
    print("=" * 70)
    print("  OK - the job can now run these steps offline." if rc == 0
          else "  INCOMPLETE - re-run on a node with internet before submitting.")
    return rc

if __name__ == "__main__":
    sys.exit(main())
