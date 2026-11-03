#!/usr/bin/env python3
"""Summarizes covariate-repro-task.sh: false-term rates from the original runs, the candidate, and covariate-refit.py.

    covariate-repro-summary.py <repro dir> <original confirm dir> <out prefix>

For each run and family, compares significant term sets. Writes <prefix>_cells.tsv. Standard library only.
"""
import csv, glob, math, os, sys
from collections import defaultdict
repro, orig, prefix = sys.argv[1:4]
csv.field_size_limit(10**9)
def r_hits(path):
    out = defaultdict(set); fams = set()
    for r in csv.DictReader(open(path), delimiter="\t", quoting=csv.QUOTE_NONE):
        if r["method"] == "covariate":
            fams.add(r["family"])
            if r["term"] != "(family)": out[r["family"]].add(r["term"])
    return out, fams
def py_hits(path):
    out = defaultdict(set)
    for r in csv.DictReader(open(path), delimiter="\t", quoting=csv.QUOTE_NONE):
        if float(r["p_adj"]) <= 0.05: out[r["family"]].add(r["term"])
    return out
def wilson(k, n, z=1.96):
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d; h = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return max(0, c-h), min(1, c+h)
cells = defaultdict(lambda: dict(n=0, rc=0, py=0, old=0, rc_py_same=0, rc_old_same=0, old_n=0, missing=0))
for f in sorted(glob.glob(f"{repro}/*/*/terms_*_*.tsv")):
    p, sc = f.split("/")[-3:-1]
    t, seed = os.path.basename(f)[6:-4].rsplit("_", 1)
    py = f"{os.path.dirname(f)}/py_{t}_{seed}.tsv"
    if not os.path.exists(py) or os.path.getsize(py) == 0:
        cells[(p, sc, t, "ALL")]["missing"] += 1; continue
    rc, fams = r_hits(f); ph = py_hits(py)
    old = None
    of = f"{orig}/{p}/{sc}/terms_{seed}.tsv"
    if t == "1" and os.path.exists(of): old = r_hits(of)[0]
    for fam in sorted(fams):
        c = cells[(p, sc, t, fam)]
        c["n"] += 1; c["rc"] += bool(rc[fam]); c["py"] += bool(ph[fam]); c["rc_py_same"] += rc[fam] == ph[fam]
        if old is not None:
            c["old_n"] += 1; c["old"] += bool(old[fam]); c["rc_old_same"] += rc[fam] == old[fam]
with open(f"{prefix}_cells.tsv", "w") as h:
    h.write("pair\tscenario\tthreshold\tfamily\truns\trate_candidate\tci_low\tci_high\trate_refit\trate_original\tcandidate_equals_refit\tcandidate_equals_original\n")
    for k in sorted(cells):
        c = cells[k]
        if not c["n"]:
            print("missing", k, c["missing"]); continue
        lo, hi = wilson(c["rc"], c["n"])
        h.write("\t".join(map(str, (*k, c["n"], f"{c['rc']/c['n']:.3f}", f"{lo:.3f}", f"{hi:.3f}", f"{c['py']/c['n']:.3f}", f"{c['old']/c['old_n']:.3f}" if c["old_n"] else "", f"{c['rc_py_same']}/{c['n']}", f"{c['rc_old_same']}/{c['old_n']}" if c["old_n"] else ""))) + "\n")
