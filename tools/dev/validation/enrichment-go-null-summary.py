#!/usr/bin/env python3
"""Summarizes summary_*.tsv and terms_*.tsv from enrichment-go-null.R.

    enrichment-go-null-summary.py <term to leave out> <rates.tsv> <by_ontology.tsv> <label>=<dir> ...

A run is one draw of dN/dS under one null, across all three ontologies. Under these nulls every
significant term is false, so the share of runs with any significant term is both the family-wise
error rate and the false discovery rate, for that run's family; it says nothing about the false
discovery rate when real associations exist. Rates are given for p_adj as written (the default,
one family across the ontologies), for BH within each ontology, and for the default with one term
left out, to show how much of a rate comes from it. Rates include a 95% Wilson interval.
"""
import csv
import glob
import math
import sys
from collections import Counter, defaultdict


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * max(0.0, c - h), 100 * min(1.0, c + h)


def read(pattern):
    rows = []
    for name in sorted(glob.glob(pattern)):
        with open(name, newline="") as handle:
            rows += list(csv.DictReader(handle, delimiter="\t"))
    return rows


def rate(k, n):
    low, high = wilson(k, n)
    return [k, f"{100 * k / n:.1f}", f"{low:.1f}", f"{high:.1f}"]


left_out, rates_path, by_ontology_path = sys.argv[1:4]
with open(rates_path, "w", newline="") as a, open(by_ontology_path, "w", newline="") as b:
    rates = csv.writer(a, delimiter="\t", lineterminator="\n")
    by_ont = csv.writer(b, delimiter="\t", lineterminator="\n")
    rates.writerow(["label", "null", "runs", "mean_raw_p05",
                    "default_any", "default_pct", "default_ci_low", "default_ci_high", "default_mean_terms",
                    "per_ontology_any", "per_ontology_pct", "per_ontology_ci_low", "per_ontology_ci_high", "per_ontology_mean_terms",
                    f"default_without_{left_out}_any", "without_pct", "without_ci_low", "without_ci_high",
                    "runs_with_left_out_term_significant", "most_frequent_terms"])
    by_ont.writerow(["label", "null", "ontology", "runs", "tested", "mean_p05", "mean_p01", "mean_p001", "min_p",
                     "runs_default_sig", "runs_per_ontology_sig"])
    for arg in sys.argv[4:]:
        label, path = arg.split("=", 1)
        summary, terms = read(path + "/summary_*.tsv"), read(path + "/terms_*.tsv")
        runs = defaultdict(list)
        for row in summary:
            runs[(row["null"], row["task"], row["run"])].append(row)
        sig = defaultdict(list)
        for row in terms:
            if float(row["p_adj"]) <= 0.05:
                sig[(row["null"], row["task"], row["run"])].append(row["term"])
        for null in sorted({key[0] for key in runs}, key=lambda x: (x != "observed", x)):
            keys = [k for k in runs if k[0] == null]
            n = len(keys)
            any_default = sum(any(int(o["n_padj05"]) > 0 for o in runs[k]) for k in keys)
            any_per = sum(any(int(o["n_bh05"]) > 0 for o in runs[k]) for k in keys)
            mean = lambda field: sum(int(o[field]) for k in keys for o in runs[k]) / n
            without = sum(any(t != left_out for t in sig[k]) for k in keys)
            common = Counter(t for k in keys for t in sig[k]).most_common(4)
            rates.writerow([label, null, n, f"{mean('n_p05'):.1f}", *rate(any_default, n), f"{mean('n_padj05'):.3f}",
                            *rate(any_per, n), f"{mean('n_bh05'):.3f}", *rate(without, n),
                            sum(left_out in sig[k] for k in keys), ";".join(f"{t}:{c}" for t, c in common)])
            for ont in ("BP", "MF", "CC"):
                rows = [o for k in keys for o in runs[k] if o["ontology"] == ont]
                f = lambda field: sum(float(o[field]) for o in rows) / len(rows)
                by_ont.writerow([label, null, ont, len(rows), rows[0]["tested"], f"{f('n_p05'):.2f}", f"{f('n_p01'):.2f}",
                                 f"{f('n_p001'):.2f}", f"{min(float(o['min_p']) for o in rows):.3g}",
                                 sum(int(o["n_padj05"]) > 0 for o in rows), sum(int(o["n_bh05"]) > 0 for o in rows)])
