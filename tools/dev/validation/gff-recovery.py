#!/usr/bin/env python3
"""Checks what the GFF comment fix changed in append_annotations output.

    gff-recovery.py <old dnds_annot.tsv> <new dnds_annot.tsv> <truth.tsv from annotate-gff.py> <query GFF>

Compares every field of every row between the two tables, then compares the terms of each
changed row with the source annotation files, which dndsR did not produce. Exits nonzero if
a field other than the query annotation changed, or a recovered term is not in the source.
InterPro IDs that the source file lacks must appear in the raw GFF lines of that gene, because
RefSeq Note attributes also carry them.
"""
import csv
import re
import sys
from collections import Counter

csv.field_size_limit(10**9)


def read(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        return reader.fieldnames, list(reader)


def terms(value):
    return set() if value in ("", "NA") else set(value.split(";"))


old_fields, old = read(sys.argv[1])
new_fields, new = read(sys.argv[2])
_, truth_rows = read(sys.argv[3])
truth = {row["gene"]: row for row in truth_rows}
ok = old_fields == new_fields and [r["query_id"] for r in old] == [r["query_id"] for r in new]
print(f"rows old={len(old)} new={len(new)} same_columns_and_order={ok}")

changed = Counter()
counts = Counter()
extras = {}
for before, after in zip(old, new):
    fields = [k for k in old_fields if before[k] != after[k]]
    changed.update(fields)
    gene = after["query_id"]
    source = truth.get(gene)
    if source is not None:
        counts["go_matches_source_new"] += terms(after["q_go"]) == terms(source["go"])
        counts["go_matches_source_old"] += terms(before["q_go"]) == terms(source["go"])
        counts["ipr_missing_new"] += bool(terms(source["ipr"]) - terms(after["q_ipr"]))
        counts["ipr_missing_old"] += bool(terms(source["ipr"]) - terms(before["q_ipr"]))
    if "q_go" not in fields and "q_ipr" not in fields:
        continue
    counts["pairs_with_changed_terms"] += 1
    if source is None:
        counts["changed_without_source"] += 1
        continue
    go_b, go_a, go_s = terms(before["q_go"]), terms(after["q_go"]), terms(source["go"])
    ipr_b, ipr_a, ipr_s = terms(before["q_ipr"]), terms(after["q_ipr"]), terms(source["ipr"])
    counts["lost_a_term"] += bool((go_b - go_a) | (ipr_b - ipr_a))
    counts["go_identical_to_source"] += go_a == go_s
    counts["go_terms_not_in_source"] += len(go_a - go_s)
    counts["ipr_contains_source"] += ipr_s <= ipr_a
    counts["ipr_extra_terms"] += len(ipr_a - ipr_s)
    if ipr_a - ipr_s:
        extras[gene] = ipr_a - ipr_s
    counts["go_terms_recovered"] += len(go_a - go_b)
    counts["ipr_terms_recovered"] += len(ipr_a - ipr_b)

locus = re.compile(r"\b(" + "|".join(map(re.escape, extras)) + r")\b") if extras else None
in_gff = {gene: set() for gene in extras}
with open(sys.argv[4], errors="replace") as handle:
    for line in handle:
        if locus and "IPR" in line and not line.startswith("#"):
            for gene in set(locus.findall(line)):
                in_gff[gene].update(re.findall(r"IPR\d{6}", line))
counts["ipr_extra_terms_not_in_gff"] = sum(len(extras[g] - in_gff[g]) for g in extras)

print("changed fields:", dict(changed))
for key in sorted(counts):
    print(f"{key}={counts[key]}")
bad = (
    not ok
    or set(changed) - {"q_attributes", "q_go", "q_ipr", "q_pfam", "q_tigrfam", "q_cog"}
    or counts["go_terms_not_in_source"]
    or counts["ipr_extra_terms_not_in_gff"]
    or counts["lost_a_term"]
    or counts["changed_without_source"]
)
sys.exit(1 if bad else 0)
