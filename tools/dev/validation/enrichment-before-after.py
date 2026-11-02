#!/usr/bin/env python3
"""Tabulates one candidate's enrichment results before and after the multiple-testing fix.

    enrichment-before-after.py <old dir> <new dir> <observed_terms.tsv> <InterPro entry list> > before_after.tsv

Each directory holds the *_q_IPR_*enrichment.tsv tables of one run. observed_terms.tsv comes
from enrichment-calibration.R, and supplies the Westfall-Young adjusted p-values.
"""
import csv
import glob
import os
import sys

csv.field_size_limit(10**9)


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


old_dir, new_dir, observed, entries = sys.argv[1:5]
names = {row["ENTRY_AC"]: row["ENTRY_NAME"] for row in read(entries)}
wy = {row["term"]: row for row in read(observed)}
out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
out.writerow(["function", "table", "term", "name", "pos", "total", "p_value", "p_adj_old", "p_adj_new",
              "sig_old", "sig_new", "wy_permute", "wy_strata", "wy_rotate"])
for new_path in sorted(glob.glob(os.path.join(new_dir, "*_q_IPR_*enrichment.tsv"))):
    base = os.path.basename(new_path)
    table = base.split("_q_IPR_")[-1].replace("enrichment.tsv", "").rstrip("_") or "all"
    fun = "term_enrichment" if table == "all" else "ipr_enrichment"
    key = "term" if fun == "term_enrichment" else "IPR"
    before = {row[key]: row for row in read(os.path.join(old_dir, base))}
    rows = read(new_path)
    rows.sort(key=lambda row: float(row["p_adj"]))
    for row in rows:
        old = before.get(row[key])
        old_adj = float(old["p_adj"]) if old else float("nan")
        extra = wy.get(row[key], {}) if fun == "term_enrichment" else {}
        out.writerow([fun, table, row[key], names.get(row[key], ""), row["pos_count"], row["total_count"], f"{float(row['p_value']):.3g}",
                      f"{old_adj:.3g}", f"{float(row['p_adj']):.3g}", int(old_adj <= 0.05), int(float(row["p_adj"]) <= 0.05),
                      *[f"{float(extra[k]):.3g}" if k in extra else "" for k in ("wy_permute", "wy_strata", "wy_rotate")]])
