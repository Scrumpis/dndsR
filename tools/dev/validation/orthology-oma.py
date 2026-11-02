#!/usr/bin/env python3
"""Compares dndsR's reciprocal best hits with OMA's pairwise orthologs, an independent resource.

    orthology-oma.py download <OMA genome 1> <OMA genome 2> <pairs.tsv>
    orthology-oma.py compare <pairs.tsv> <dnds.tsv> <query.faa.gz> <query.gff> <subject.faa.gz> <subject.gff> <prefix>
    orthology-oma.py by_ds <prefix>_pairs.tsv <out.tsv>

OMA's identifiers differ from RefSeq's, so a protein is matched by the MD5 of its sequence, then by
its accession. Only genes that map on both sides of an OMA pair can be scored, and a query gene can
have several OMA orthologs. Each dndsR pair is then one of:

    agree     its subject is among the OMA orthologs of its query
    disagree  its query has OMA orthologs that map, and its subject is not one of them
    unscored  its query has none

Writes <prefix>_pairs.tsv (per dndsR pair) and <prefix>_summary.tsv. by_ds tabulates scored pairs by dS,
and for pairs above 3 times the median dS. Standard library only.
"""
import csv
import gzip
import hashlib
import json
import math
import re
import statistics
import sys
import time
import urllib.request
from collections import defaultdict

csv.field_size_limit(10**9)


def download(a, b, out):
    page, total = 1, 0
    with open(out, "w") as handle:
        handle.write("q_id\tq_md5\ts_id\ts_md5\trel_type\tscore\tdistance\n")
        while True:
            url = f"https://omabrowser.org/api/pairs/{a}/{b}/?per_page=1000&page={page}"
            for attempt in range(6):
                try:
                    # OMA's firewall refuses Python's default client name, so this one says who is asking.
                    request = urllib.request.Request(url, headers={"User-Agent": "dndsR-validation (https://github.com/Scrumpis/dndsR)"})
                    with urllib.request.urlopen(request, timeout=180) as r:
                        data = json.load(r)
                    break
                except Exception as err:
                    print(f"retry page {page}: {err}", file=sys.stderr)
                    time.sleep(30 * (attempt + 1))
            else:
                sys.exit(f"could not read page {page}")
            if not data:
                break
            for row in data:
                x, y = row["entry_1"], row["entry_2"]
                handle.write(f"{x['canonicalid']}\t{x['sequence_md5']}\t{y['canonicalid']}\t{y['sequence_md5']}"
                             f"\t{row.get('rel_type')}\t{row.get('score')}\t{row.get('distance')}\n")
            total += len(data)
            print(f"page {page}: {total} pairs", file=sys.stderr)
            page += 1
            time.sleep(1)


def protein_to_gene(faa, gff):
    """{protein accession: gene}, {sequence MD5: set of genes}. Genes named as in dndsR's tables."""
    accession = {}
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] == "CDS":
                attrs = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
                gene = attrs.get("locus_tag") or attrs.get("gene")
                if "protein_id" in attrs and gene:
                    accession[attrs["protein_id"]] = gene
    by_md5 = defaultdict(set)
    name, parts = None, []

    def flush():
        if name in accession:
            by_md5[hashlib.md5("".join(parts).encode()).hexdigest()].add(accession[name])

    with gzip.open(faa, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip().rstrip("*"))
        flush()
    return accession, by_md5


def lookup(identifier, md5, accession, by_md5):
    genes = set(by_md5.get(md5, ()))
    if identifier in accession:
        genes.add(accession[identifier])
    return genes


def number(x):
    try:
        v = float(x)
    except ValueError:
        return math.nan
    return v


def compare(pairs, dnds, qfaa, qgff, sfaa, sgff, prefix):
    qacc, qmd5 = protein_to_gene(qfaa, qgff)
    sacc, smd5 = protein_to_gene(sfaa, sgff)
    orthologs = defaultdict(set)
    total = mapped = 0
    with open(pairs) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            total += 1
            q = lookup(row["q_id"], row["q_md5"], qacc, qmd5)
            s = lookup(row["s_id"], row["s_md5"], sacc, smd5)
            if q and s:
                mapped += 1
                for g in q:
                    orthologs[g] |= s
    rows = []
    with open(dnds) as handle:
        for row in csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE):
            q, s = row["query_id"], row["subject_id"]
            status = "unscored" if q not in orthologs else ("agree" if s in orthologs[q] else "disagree")
            rows.append((q, s, status, number(row["dN"]), number(row["dS"]), number(row["dNdS"])))
    with open(f"{prefix}_pairs.tsv", "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["query_id", "subject_id", "status", "dN", "dS", "dNdS"])
        w.writerows(rows)
    with open(f"{prefix}_summary.tsv", "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["status", "pairs", "percent", "median_dS", "dS_above_1_pct", "median_dNdS", "pos_gt1", "pos_pct", "oma_pairs", "oma_pairs_mapped", "oma_query_genes"])
        for status in ("agree", "disagree", "unscored"):
            sel = [r for r in rows if r[2] == status]
            ds = [r[4] for r in sel if r[4] == r[4]]
            ok = [r[5] for r in sel if r[5] == r[5] and r[5] < 10]
            w.writerow([status, len(sel), f"{100 * len(sel) / len(rows):.1f}",
                        f"{statistics.median(ds):.3f}" if ds else "NA",
                        f"{100 * sum(x > 1 for x in ds) / max(1, len(ds)):.1f}",
                        f"{statistics.median(ok):.3f}" if ok else "NA",
                        sum(x > 1 for x in ok), f"{100 * sum(x > 1 for x in ok) / max(1, len(ok)):.2f}",
                        total, mapped, len(orthologs)])


def by_ds(pairs, out):
    with open(pairs) as handle:
        rows = [(r["status"], number(r["dS"])) for r in csv.DictReader(handle, delimiter="\t")]
    known = [d for _, d in rows if d == d]
    cut = 3 * statistics.median(known)
    classes = [("dS below 0.02", lambda d: d < 0.02), ("0.02 to 0.1", lambda d: 0.02 <= d < 0.1), ("0.1 to 0.3", lambda d: 0.1 <= d < 0.3),
               ("0.3 to 0.8", lambda d: 0.3 <= d < 0.8), ("0.8 to 1.5", lambda d: 0.8 <= d < 1.5), ("1.5 and above", lambda d: d >= 1.5),
               (f"above 3 x median ({cut:.3f})", lambda d: d > cut), (f"at most 3 x median ({cut:.3f})", lambda d: d <= cut)]
    with open(out, "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["class", "agree", "disagree", "unscored", "disagree_pct_of_scored"])
        for label, test in classes:
            n = {k: sum(s == k and d == d and test(d) for s, d in rows) for k in ("agree", "disagree", "unscored")}
            scored = n["agree"] + n["disagree"]
            w.writerow([label, n["agree"], n["disagree"], n["unscored"], f"{100 * n['disagree'] / scored:.1f}" if scored else "NA"])


if sys.argv[1] == "download":
    download(*sys.argv[2:5])
elif sys.argv[1] == "by_ds":
    by_ds(*sys.argv[2:4])
else:
    compare(*sys.argv[2:9])
