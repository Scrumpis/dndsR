#!/usr/bin/env python3
"""Summarizes real_null_<seed>.tsv files from enrichment-real-null.R.

    enrichment-null-summary.py <label>=<dir> ...

A run is one draw for one function under one null. Its significant terms are summed across
the tables it wrote, and every one of them is a false positive. Runs with a recorded error
are left out and counted. Rates include a 95% Wilson interval.
"""
import csv
import glob
import math
import sys
from collections import defaultdict


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


print("label\tfun\tnull\truns\terrors\tany_sig\tany_sig_pct\tci_low\tci_high\tmean_sig\tmean_tested_per_table")
for arg in sys.argv[1:]:
    label, path = arg.split("=", 1)
    runs = defaultdict(lambda: {"sig": 0, "tested": [], "error": False})
    for name in sorted(glob.glob(path + "/real_null_*.tsv")):
        with open(name, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                run = runs[(row["fun"], row["null"], name, row["run"])]
                run["sig"] += int(row["n_sig"])
                run["error"] = run["error"] or row["file"] == "(error)"
                if row["n_tested"] not in ("", "NA"):
                    run["tested"].append(int(row["n_tested"]))
    groups = defaultdict(list)
    for (fun, null, _, _), run in runs.items():
        groups[(fun, null)].append(run)
    for (fun, null), group in sorted(groups.items()):
        done = [r for r in group if not r["error"]]
        n = len(done)
        k = sum(r["sig"] > 0 for r in done)
        low, high = wilson(k, n)
        tested = [t for r in done for t in r["tested"]]
        mean_tested = sum(tested) / len(tested) if tested else float("nan")
        print(f"{label}\t{fun}\t{null}\t{n}\t{len(group) - n}\t{k}\t{100 * k / n:.1f}\t{100 * low:.1f}\t{100 * high:.1f}"
              f"\t{sum(r['sig'] for r in done) / n:.2f}\t{mean_tested:.0f}")
