#!/usr/bin/env python3
"""Asks whether one GO term's excess under a stratified null falls as the strata get finer.

    enrichment-go-confounding.py <dnds_annot.tsv> <ancestors.tsv> <term> <p threshold> <draws> <seed>

Reassigns dN/dS within strata of q_len x dS x GO richness, at several numbers of quantile bins,
and recalculates the term's classic one-sided Fisher p-value from counts. The term's conditional
weight01 p-value is not used, so these rates describe confounding, not the default's error rate.
A rate falling with finer strata means that residual length or dS differences produce the excess.
Standard library only; descriptive, and no input to dndsR.
"""
import csv
import math
import random
import statistics
import sys

csv.field_size_limit(10**9)


def lchoose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def upper_tail(a, m, k, n):
    high = min(m, k)
    return min(1.0, sum(math.exp(lchoose(k, x) + lchoose(n - k, m - x) - lchoose(n, m)) for x in range(a, high + 1)))


def bins(values, k):
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0] * len(values)
    for rank, i in enumerate(order):
        out[i] = rank * k // len(values)
    return out


annot, ancestors, term, cut, draws, seed = sys.argv[1:7]
cut, draws = float(cut), int(draws)
random.seed(int(seed))
ontology, parents = {}, {}
with open(ancestors, newline="") as h:
    for r in csv.DictReader(h, delimiter="\t", quoting=csv.QUOTE_NONE):
        parents.setdefault(r["term"], set()).add(r["ancestor"])
        ontology[r["term"]] = r["ontology"]
ont = ontology[term]
rows = []
with open(annot, newline="") as h:
    for r in csv.DictReader(h, delimiter="\t", quoting=csv.QUOTE_NONE):
        try:
            ratio, ds, length = float(r["dNdS"]), float(r["dS"]), float(r["q_len"])
        except ValueError:
            continue
        if math.isfinite(ratio) and ratio < 10:
            go = [t for t in r["q_go"].split(";") if t and t != "NA"]
            rows.append((ratio, ds, length, go))
n = len(rows)
ratio = [r[0] for r in rows]
ds = [r[1] for r in rows]
length = [r[2] for r in rows]
counts = [len(r[3]) for r in rows]
middle = statistics.median(c for c in counts if c)
rich = [0 if c == 0 else 1 if c <= middle else 2 for c in counts]
in_ont = [any(ontology.get(t) == ont for t in r[3]) for r in rows]
member = [any(term in parents.get(t, ()) for t in r[3] if ontology.get(t) == ont) for r in rows]
universe = sum(in_ont)
m = sum(member)


def pvalue(values):
    k = sum(1 for v, u in zip(values, in_ont) if u and v > 1)
    a = sum(1 for v, t in zip(values, member) if t and v > 1)
    return a, k, upper_tail(a, m, k, universe)


a, k, p = pvalue(ratio)
inside = [i for i in range(n) if member[i]]
others = [i for i in range(n) if in_ont[i] and not member[i]]
print(f"term\t{term}\tontology\t{ont}\tuniverse\t{universe}\tmembers\t{m}")
print(f"observed\tpositive_in_term\t{a}\tpositive_in_universe\t{k}\tclassic_p\t{p:.3g}")
print("median\tq_len_term\t%g\tq_len_rest\t%g\tdS_term\t%.4f\tdS_rest\t%.4f" % (
    statistics.median(length[i] for i in inside), statistics.median(length[i] for i in others),
    statistics.median(ds[i] for i in inside), statistics.median(ds[i] for i in others)))
print("bins\tstrata\tdraws\tmean_positive_in_term\tpct_p_lt_0.05\tpct_p_lt_cut\tdraws_as_extreme_as_observed")
for b in (5, 10, 20, 40):
    lb, db = bins(length, b), bins(ds, b)
    groups = {}
    for i in range(n):
        groups.setdefault((lb[i], db[i], rich[i]), []).append(i)
    hits05 = hits = extreme = total = 0
    for _ in range(draws):
        x = list(ratio)
        for s in groups.values():
            v = [ratio[i] for i in s]
            random.shuffle(v)
            for i, w in zip(s, v):
                x[i] = w
        got, _, q = pvalue(x)
        total += got
        hits05 += q < 0.05
        hits += q < cut
        extreme += q <= p
    print(f"{b}\t{len(groups)}\t{draws}\t{total / draws:.1f}\t{100 * hits05 / draws:.1f}\t{100 * hits / draws:.1f}\t{extreme}")
