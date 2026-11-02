#!/usr/bin/env python3
"""Recalculates covariate_enrichment from a table with no dndsR or R code, including the logistic fit.

    covariate-refit.py <N_dnds_annot.tsv> <go_ancestors.tsv> <out.tsv> [threshold] [R table prefix]

Writes one row per tested term: type, ontology, term, total, pos, expected, p, p_adj. With an R table prefix,
such as dir/N, it compares <prefix>_q_GO_covariate_enrichment.tsv and _IPR_, and prints the largest differences.

Written from the method as documented in docs/methods.md, not from R/covariate_enrichment.R's code:
- natural cubic splines from the truncated-power basis in Hastie, Tibshirani, and Friedman, Elements of
  Statistical Learning, eq. 5.4-5.5. Fitted values do not depend on which basis spans the same space.
  Boundary knots at the range and interior knots at the 1/3 and 2/3 quantiles, 3 degrees of freedom.
- logistic regression by Newton steps on a QR factorization
- Poisson-binomial tail by direct convolution
- Benjamini-Hochberg by the step-up definition
Needs numpy. go_ancestors.tsv has columns ontology, term, ancestor, as enrichment-go-ancestors.R writes it.
"""
import csv
import sys
from collections import defaultdict

import numpy as np

csv.field_size_limit(10**9)
ROOTS = {"BP": "GO:0008150", "MF": "GO:0003674", "CC": "GO:0005575"}


def natural_basis(x):
    lo, hi = x.min(), x.max()
    k1, k2 = np.quantile(x, [1 / 3, 2 / 3])  # the default definition is R's type 7
    knots = [lo, k1, k2, hi]

    def d(k):
        return (np.maximum(x - k, 0) ** 3 - np.maximum(x - hi, 0) ** 3) / (hi - k)

    return np.column_stack([x, d(lo) - d(k2), d(k1) - d(k2)]), knots


def fit(y, length, ds):
    a, _ = natural_basis(np.log(length))
    b, _ = natural_basis(np.log(ds))
    X = np.column_stack([np.ones(len(y)), a, b])
    X[:, 1:] = (X[:, 1:] - X[:, 1:].mean(0)) / X[:, 1:].std(0)
    Q, R = np.linalg.qr(X)
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean() / (1 - y.mean()))
    last = np.inf
    for _ in range(100):
        eta = X @ beta
        mu = np.clip(1 / (1 + np.exp(-eta)), 1e-15, 1 - 1e-15)
        dev = -2 * np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu))
        if abs(dev - last) < 1e-11 * (abs(dev) + 0.1):
            break
        last = dev
        w = mu * (1 - mu)
        z = eta + (y - mu) / w
        sw = np.sqrt(w)
        beta = np.linalg.lstsq(X * sw[:, None], z * sw, rcond=None)[0]
    return 1 / (1 + np.exp(-(X @ beta)))


def tail(p, k):
    if k <= 0:
        return 1.0
    if k > len(p):
        return 0.0
    f = np.zeros(k + 1)
    f[0] = 1.0
    for q in p:
        absorbed = f[k] + f[k - 1] * q
        f[1:k] = f[1:k] * (1 - q) + f[:k - 1] * q
        f[0] *= 1 - q
        f[k] = absorbed
    return min(1.0, f[k])


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p, kind="stable")
    adj = p[order] * len(p) / np.arange(1, len(p) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(adj, 1)
    return out


def universe(gene_terms, dnds, length, ds, threshold, min_total=10):
    has = [i for i, t in enumerate(gene_terms) if t]
    if not has:
        return []
    d = dnds[has]
    cut = threshold if not isinstance(threshold, str) else np.quantile(d, float(threshold.lstrip("q")))
    y = (d > cut).astype(float)
    rate = y.mean()
    prob = np.full(len(has), rate) if rate in (0.0, 1.0) else np.maximum(fit(y, length[has], ds[has]), 0.1 * rate)
    members = defaultdict(list)
    for j, i in enumerate(has):
        for t in gene_terms[i]:
            members[t].append(j)
    rows = []
    for t, m in members.items():
        if len(m) >= min_total:
            k = int(y[m].sum())
            rows.append((t, len(m), k, float(prob[m].sum()), tail(prob[m], k)))
    return rows


def main():
    table, ancestors, out = sys.argv[1:4]
    threshold = sys.argv[4] if len(sys.argv) > 4 else "1"
    threshold = threshold if threshold.startswith("q") else float(threshold)
    prefix = sys.argv[5] if len(sys.argv) > 5 else None
    with open(table, newline="") as h:
        rows = list(csv.DictReader(h, delimiter="\t", quoting=csv.QUOTE_NONE))

    def num(x):
        try:
            return float(x)
        except ValueError:
            return float("nan")

    keep = [r for r in rows if num(r["dNdS"]) < 10 and num(r["dS"]) > 0 and num(r["q_len"]) > 0]
    dnds = np.array([num(r["dNdS"]) for r in keep])
    ds = np.array([num(r["dS"]) for r in keep])
    length = np.array([num(r["q_len"]) for r in keep])

    def terms(col):
        return [{t for t in r[col].replace(",", ";").replace("|", ";").split(";") if t and t != "NA"} for r in keep]

    results = []
    ipr = universe(terms("q_ipr"), dnds, length, ds, threshold)
    results += [("IPR", "", *x) for x in ipr]
    close = defaultdict(lambda: defaultdict(set))
    with open(ancestors, newline="") as h:
        for a in csv.DictReader(h, delimiter="\t"):
            close[a["ontology"]][a["term"]].add(a["ancestor"])
    go_terms = terms("q_go")
    go_rows = []
    for o in ("BP", "MF", "CC"):
        known = close[o]
        closed = [set().union(*[(known[t] | {t}) - {"all"} for t in g if t in known], set()) | ({ROOTS[o]} if ROOTS[o] in g else set()) for g in go_terms]
        go_rows += [("GO", o, *x) for x in universe(closed, dnds, length, ds, threshold)]
    results += go_rows
    adj = {}
    for family in ("IPR", "GO"):
        idx = [i for i, r in enumerate(results) if r[0] == family]
        if idx:
            for i, q in zip(idx, bh([results[i][6] for i in idx])):
                adj[i] = float(q)
    with open(out, "w") as h:
        h.write("family\tontology\tterm\ttotal\tpos\texpected\tp\tp_adj\n")
        for i, r in enumerate(results):
            h.write("\t".join(map(str, r)) + f"\t{adj[i]}\n")
    if prefix:
        mine = {(r[0], r[1], r[2]): (r[3], r[4], r[5], r[6], adj[i]) for i, r in enumerate(results)}
        for family in ("GO", "IPR"):
            with open(f"{prefix}_q_{family}_covariate_enrichment.tsv", newline="") as h:
                theirs = list(csv.DictReader(h, delimiter="\t", quoting=csv.QUOTE_NONE))
            matched = 0; worst = [0.0, 0.0, 0.0, 0.0]; counts = 0; hits = 0
            for t in theirs:
                key = (family, t["ontology"] if family == "GO" else "", t["term"])
                if key not in mine:
                    counts += 1
                    continue
                total, pos, exp, p, q = mine[key]
                matched += 1
                counts += (total != int(t["total_count"])) + (pos != int(t["pos_count"]))
                worst[0] = max(worst[0], abs(exp - float(t["expected"])) / max(exp, 1e-12))
                lp = lambda v: np.log10(max(v, 1e-300))
                worst[1] = max(worst[1], abs(lp(p) - lp(float(t["p_value"]))))
                worst[2] = max(worst[2], abs(p - float(t["p_value"])))
                worst[3] = max(worst[3], abs(q - float(t["p_adj"])))
                hits += (q <= 0.05) != (float(t["p_adj"]) <= 0.05)
            print(f"{family}\tcompared={matched}\tof_python={sum(1 for k in mine if k[0] == family)}\tcount_differences={counts}\tcalls_differ={hits}"
                  f"\tmax_rel_expected={worst[0]:.3g}\tmax_log10p={worst[1]:.3g}\tmax_abs_p={worst[2]:.3g}\tmax_abs_padj={worst[3]:.3g}")


main()
