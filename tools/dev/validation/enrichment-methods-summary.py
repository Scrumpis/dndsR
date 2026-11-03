#!/usr/bin/env python3
"""Summarizes enrichment-methods.R over null and spiked runs.

    enrichment-methods-summary.py <runs dir> <dnds_annot.tsv> <go_ancestors.tsv> <out prefix>

<runs dir>/<scenario>/terms_<seed>.tsv and, for simulations, sim_<seed>.tsv.gz.

Null scenarios have no true association, so any significant term in a run is false: the share of
runs with one is the family-wise error rate and the false discovery rate under the complete null.

Spiked scenarios raise omega in a known set of genes, in 5 chosen InterPro terms. A significant term
counts as:
    spiked   it is one of the 5
    related  its genes overlap the spiked genes more than chance (hypergeometric p < 0.001), as a term
             sharing genes with a spiked one does
    false    otherwise
The false discovery proportion of a run is false / (all significant), 0 for a run with none, and the
false discovery rate is its mean. Power is the share of spiked InterPro terms found. Standard library only.
"""
import csv
import glob
import gzip
import math
import os
import sys
from collections import defaultdict
from math import comb

csv.field_size_limit(10**9)
runs_dir, annot, ancestors, prefix = sys.argv[1:5]


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * max(0.0, c - h), 100 * min(1.0, c + h))


# Term membership, read only if spiked scenarios need it.
members = None


def load_members():
    global members, all_genes
    up = defaultdict(set)
    with open(ancestors) as h:
        for r in csv.DictReader(h, delimiter="\t"):
            up[r["term"]].add(r["ancestor"])
    members = defaultdict(set)
    all_genes = set()
    with open(annot, newline="") as h:
        for r in csv.DictReader(h, delimiter="\t", quoting=csv.QUOTE_NONE):
            g = r["query_id"]
            all_genes.add(g)
            for t in filter(None, r["q_ipr"].replace(",", ";").split(";")):
                if t != "NA":
                    members[t].add(g)
            for t in filter(None, r["q_go"].replace(",", ";").split(";")):
                if t != "NA":
                    for a in up.get(t, (t,)):
                        members[a].add(g)


def overlap_p(k, term, spiked):
    """P(overlap >= k) by chance, among the annotated genes."""
    universe = len(all_genes)
    if k == 0:
        return 1.0
    denom = comb(universe, spiked)
    return min(1.0, sum(comb(term, x) * comb(universe - term, spiked - x) for x in range(k, min(term, spiked) + 1)) / denom)


def read_sim(path):
    spiked_terms, spiked_genes = [], set()
    with gzip.open(path, "rt") as h:
        first = h.readline()
        spiked_terms = [x for x in first.strip().split("spiked_terms=")[-1].split(",") if x.startswith("IPR")]
        for r in csv.DictReader(h, delimiter="\t"):
            if r["spiked"] == "TRUE":
                spiked_genes.add(r["query_id"])
    return spiked_terms, spiked_genes


rows, power_rows = [], []
for folder in sorted(glob.glob(os.path.join(runs_dir, "*"))):
    scenario = os.path.basename(folder)
    spike = scenario.startswith("spike_")
    if spike and members is None:
        load_members()
    tally = defaultdict(lambda: dict(runs=0, any=0, terms=0, fdp=0.0, found=0, spiked=0, raw=0))
    frequent = defaultdict(lambda: defaultdict(int))
    for path in sorted(glob.glob(os.path.join(folder, "terms_*.tsv"))):
        seed = os.path.basename(path)[6:-4]
        with open(path, newline="") as h:
            table = list(csv.DictReader(h, delimiter="\t"))
        truth = read_sim(os.path.join(folder, f"sim_{seed}.tsv.gz")) if spike else None
        groups = defaultdict(list)
        methods = {(r["family"], r["method"]) for r in table}
        for r in table:
            if r["term"] == "(family)":
                tally[(r["family"], r["method"])]["raw"] += int(r["n_raw05"])
            else:
                groups[(r["family"], r["method"])].append(r)
        for key in methods:
            t = tally[key]
            hits = groups[key]
            t["runs"] += 1
            t["any"] += bool(hits)
            for r in hits:
                frequent[key][r["term"]] += 1
            if spike:
                terms = [r["term"] for r in hits]
                def is_false(term):
                    if term in truth[0]:
                        return False
                    m = members.get(term, set())
                    return not m or overlap_p(len(m & truth[1]), len(m), len(truth[1])) >= 0.001
                n_false = sum(is_false(x) for x in terms)
                t["terms"] += n_false
                t["fdp"] += n_false / len(terms) if terms else 0.0
                if key[0] == "IPR":
                    t["found"] += sum(x in terms for x in truth[0])
                    t["spiked"] += len(truth[0])
            else:
                t["terms"] += len(hits)
    for key, t in sorted(tally.items()):
        n = t["runs"]
        low, high = wilson(t["any"], n)
        top = sorted(frequent[key].items(), key=lambda x: -x[1])[:3]
        rows.append([scenario, key[0], key[1], n, f"{t['raw'] / max(n, 1):.1f}", t["any"], f"{100 * t['any'] / max(n, 1):.1f}", f"{low:.1f}", f"{high:.1f}",
                     f"{t['terms'] / max(n, 1):.3f}", f"{100 * t['fdp'] / max(n, 1):.1f}" if spike else "",
                     f"{100 * t['found'] / t['spiked']:.1f}" if t["spiked"] else "", ";".join(f"{a}:{b}" for a, b in top)])

with open(prefix + ".tsv", "w", newline="") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    w.writerow(["scenario", "family", "method", "runs", "mean_raw_p05", "runs_with_term", "pct", "ci_low", "ci_high",
                "mean_false_terms", "fdr_pct", "power_pct", "most_frequent_terms"])
    w.writerows(rows)
