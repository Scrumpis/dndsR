#!/usr/bin/env python3
"""Summarizes estimator-benchmark.R: bias, spread, failure, calls above 1, and bootstrap coverage.

    estimator-benchmark-summary.py <dir of cell_*.tsv> <out prefix>

<prefix>_cells.tsv has a row per estimator and cell. For each it gives the median of estimate / truth for the
ratio and dS, the geometric spread (median absolute log error of the ratio), the share undefined,
and the share above 1. Comeron also has the percentile bootstrap's coverage of the true ratio, and
how often its lower bound is above 1. Standard library only.
"""
import csv
import glob
import math
import statistics
import sys

folder, prefix = sys.argv[1:3]
methods = ["Comeron", "YN00", "NG86", "ML", "MLWL", "LPB"]


def f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return math.nan
    return v


with open(prefix + "_cells.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["regime", "codons", "true_dS", "true_omega", "true_dS_actual", "method", "reps", "undefined_pct", "ratio_median_over_truth",
                "ratio_median_abs_log_error", "dS_median_over_truth", "dN_median_over_truth", "above1_pct", "boot_n", "boot_coverage_pct", "boot_lower_above1_pct", "boot_median_width_ratio"])
    for path in sorted(glob.glob(folder + "/cell_*.tsv")):
        with open(path, newline="") as h:
            rows = list(csv.DictReader(h, delimiter="\t"))
        if not rows:
            continue
        r0 = rows[0]
        true_dn, true_ds = f(r0["true_dN"]), f(r0["true_dS"])
        truth = true_dn / true_ds
        for m in methods:
            dn = [f(r[m + "_dN"]) for r in rows]
            ds = [f(r[m + "_dS"]) for r in rows]
            good = [(a, b) for a, b in zip(dn, ds) if a == a and b == b and b > 0 and b < 90 and a < 90]
            ratios = [a / b for a, b in good]
            pos = [x for x in ratios if x > 0]
            bn = cov = low = 0
            widths = []
            if m == "Comeron":
                for r in rows:
                    lo, hi = f(r["boot_lo"]), f(r["boot_hi"])
                    if lo == lo and hi == hi:
                        bn += 1
                        cov += lo <= truth <= hi
                        low += lo > 1
                        if lo > 0:
                            widths.append(hi / lo)
            w.writerow([r0["regime"], r0["codons"], r0["ds"], r0["omega"], f"{true_ds:.4f}", m, len(rows),
                        f"{100 * (1 - len(good) / len(rows)):.1f}",
                        f"{statistics.median(ratios) / truth:.3f}" if ratios else "NA",
                        f"{statistics.median(abs(math.log(x / truth)) for x in pos):.3f}" if pos else "NA",
                        f"{statistics.median(b for _, b in good) / true_ds:.3f}" if good else "NA",
                        f"{statistics.median(a for a, _ in good) / true_dn:.3f}" if good else "NA",
                        f"{100 * sum(x > 1 for x in ratios) / max(1, len(ratios)):.1f}",
                        bn if bn else "", f"{100 * cov / bn:.1f}" if bn else "", f"{100 * low / bn:.1f}" if bn else "",
                        f"{statistics.median(widths):.2f}" if widths else ""])
