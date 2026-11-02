#!/usr/bin/env python3
"""Recalculate dndsR's Fisher p-values and BH from a per-gene table, without R or dndsR.

    enrichment-verify.py term <dnds_annot.tsv> <enrichment.tsv> [min_total] [min_pos]
    enrichment-verify.py ipr <dnds_annot.tsv> <entry list> <dir with tables> <prefix> [min_total] [min_pos] [max_prop]

Prints each table's family size, rows, and largest differences from the independent values,
and exits nonzero if a table disagrees. Standard library only.
"""
import csv
import glob
import math
import os
import sys
from collections import Counter, defaultdict

csv.field_size_limit(10**9)


def lchoose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def upper_tail(a, term, positive, universe):
    """P(X >= a) for the hypergeometric distribution, equal to one-sided Fisher's test."""
    low, high = max(a, term + positive - universe, 0), min(term, positive)
    denom = lchoose(universe, positive)
    total = sum(math.exp(lchoose(term, x) + lchoose(universe - term, positive - x) - denom) for x in range(low, high + 1))
    return min(1.0, total)


def bh(pvalues):
    """Benjamini-Hochberg, from the definition."""
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    out, best = [0.0] * n, 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        best = min(best, pvalues[i] * n / rank)
        out[i] = best
    return out


def read_genes(path, column="q_ipr", max_dnds=10.0):
    genes = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE):
            try:
                ratio = float(row["dNdS"])
            except ValueError:
                continue
            if math.isfinite(ratio) and ratio < max_dnds:
                raw = row.get(column, "")
                terms = set(t for t in raw.split(";") if t and t != "NA")
                genes.append((ratio, terms))
    return genes


def family(genes, min_total, max_prop=None):
    """Annotated genes only; a term's outcome plays no part in choosing what is tested."""
    kept = [g for g in genes if g[1]]
    universe, positive = len(kept), sum(1 for g in kept if g[0] > 1)
    total, pos = Counter(), Counter()
    for ratio, terms in kept:
        total.update(terms)
        if ratio > 1:
            pos.update(terms)
    tested = [t for t in sorted(total) if total[t] >= min_total and (max_prop is None or total[t] / universe <= max_prop)]
    return {t: (pos[t], total[t], upper_tail(pos[t], total[t], positive, universe)) for t in tested}


def compare(label, fam, table, key, min_pos):
    p = {t: v[2] for t, v in fam.items()}
    names = list(p)
    adj = dict(zip(names, bh([p[t] for t in names])))
    with open(table, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))
    want = {t for t, v in fam.items() if v[0] >= min_pos}
    got = {r[key] for r in rows}
    dp = max((abs(float(r["p_value"]) - p[r[key]]) / max(p[r[key]], 1e-300) for r in rows), default=0.0)
    da = max((abs(float(r["p_adj"]) - adj[r[key]]) / max(adj[r[key]], 1e-300) for r in rows), default=0.0)
    same_rows = want == got
    n_ok = all(int(r["n_tested"]) == len(fam) for r in rows)
    print(f"{label}\tfamily={len(fam)}\trows={len(rows)}\tsignificant={sum(float(r['p_adj']) <= 0.05 for r in rows)}"
          f"\tsame_rows={same_rows}\tn_tested_ok={n_ok}\tmax_rel_dp={dp:.2e}\tmax_rel_dadj={da:.2e}")
    return same_rows and n_ok and dp < 1e-8 and da < 1e-8


def main():
    mode = sys.argv[1]
    ok = True
    if mode == "term":
        annot, table = sys.argv[2:4]
        min_total, min_pos = [int(x) for x in (sys.argv[4:6] or ["2", "2"])]
        ok = compare(os.path.basename(table), family(read_genes(annot), min_total), table, "term", min_pos)
    else:
        annot, entries, folder, prefix = sys.argv[2:6]
        rest = sys.argv[6:]
        min_total, min_pos = (int(rest[0]), int(rest[1])) if rest else (2, 2)
        max_prop = float(rest[2]) if len(rest) > 2 else 0.2
        kind = {}
        with open(entries, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE):
                kind[row["ENTRY_AC"]] = row["ENTRY_TYPE"]
        genes = read_genes(annot)
        for type_name in sorted(set(kind.values())):
            typed = [(r, {t for t in terms if kind.get(t) == type_name}) for r, terms in genes]
            fam = family(typed, min_total, max_prop)
            table = os.path.join(folder, f"{prefix}_IPR_{type_name}_enrichment.tsv")
            if os.path.exists(table):
                ok = compare(type_name, fam, table, "IPR", min_pos) and ok
            elif any(v[0] >= min_pos for v in fam.values()):
                print(f"{type_name}\tMISSING TABLE")
                ok = False
    sys.exit(0 if ok else 1)


main()
