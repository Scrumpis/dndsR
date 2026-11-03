#!/usr/bin/env python3
"""Recalculates a covariate_enrichment table without dndsR: term counts, Fisher's p-values, BH, and tail bounds.

    covariate-verify.py <dnds_annot.tsv> <table.tsv> <go_ancestors.tsv or -> <side prefix q|s>

The spline fit is not recalculated. It is checked by what must hold: expected counts lie between the floor
and the term size, and adjusted p-values follow from the p-values. Standard library only.
"""
import csv
import math
import sys
from collections import defaultdict
from math import lgamma

csv.field_size_limit(10**9)
annot, table, ancestors, side = sys.argv[1:5]


def log_choose(n, k):
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def hyper_upper(k, K, N, n):
    if k <= 0:
        return 1.0
    return min(1.0, sum(math.exp(log_choose(K, x) + log_choose(N - K, n - x) - log_choose(N, n)) for x in range(k, min(K, n) + 1)))


with open(table, newline="") as h:
    rows = list(csv.DictReader(h, delimiter="\t", quoting=csv.QUOTE_NONE))
go = rows[0]["term_type"] == "GO"
up = defaultdict(lambda: defaultdict(set))
if go:
    with open(ancestors) as h:
        for r in csv.DictReader(h, delimiter="\t"):
            up[r["ontology"]][r["term"]].add(r["ancestor"])
col = f"{side}_{'go' if go else 'ipr'}"
genes = []
with open(annot, newline="") as h:
    for r in csv.DictReader(h, delimiter="\t", quoting=csv.QUOTE_NONE):
        try:
            ratio, ds, length = float(r["dNdS"]), float(r["dS"]), float(r[f"{side}_len"])
        except ValueError:
            continue
        if not (math.isfinite(ratio) and ratio < 10 and math.isfinite(ds) and ds > 0 and length > 0):
            continue
        terms = {t for t in r[col].replace(",", ";").replace("|", ";").split(";") if t and t != "NA"}
        genes.append((ratio > 1, terms))


def universe(ontology):
    members = defaultdict(int), defaultdict(int)
    N = K = 0
    for positive, terms in genes:
        if go:
            terms = set().union(*[up[ontology].get(t, set()) for t in terms]) if terms else set()
        if terms:
            N += 1
            K += positive
            for t in terms:
                members[0][t] += 1
                members[1][t] += positive
    return N, K, members


cache = {}
bad = 0
for r in rows:
    key = r["ontology"] or "IPR"
    cache.setdefault(key, universe(key))
    N, K, (total, pos) = cache[key]
    t = r["term"]
    expected = (N, K, total[t], pos[t])
    got = tuple(int(float(r[x])) for x in ("universe", "pos_total", "total_count", "pos_count"))
    fisher = hyper_upper(pos[t], K, N, total[t])
    floor = 0.1 * K / N * total[t]
    if got != expected or abs(fisher - float(r["fisher_p"])) > 1e-6 * max(fisher, 1e-12) + 1e-14 or not floor * (1 - 1e-9) <= float(r["expected"]) <= total[t]:
        bad += 1
        if bad <= 5:
            print("DIFFERS", t, key, got, expected, fisher, r["fisher_p"], r["expected"])
n_tested = int(rows[0]["n_tested"])
if len(rows) == n_tested:
    order = sorted(range(len(rows)), key=lambda i: float(rows[i]["p_value"]))
    running = 1.0
    for rank in range(len(order) - 1, -1, -1):
        i = order[rank]
        running = min(running, float(rows[i]["p_value"]) * n_tested / (rank + 1))
        if abs(running - float(rows[i]["p_adj"])) > 1e-6 * running + 1e-14:
            bad += 1
print(f"{table}: {len(rows)} rows of {n_tested} tested; differences {bad}; BH {'checked' if len(rows) == n_tested else 'not checked: rows filtered'}")
sys.exit(bad > 0)
