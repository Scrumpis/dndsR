#!/usr/bin/env python3
"""Adds UniProt GO and InterPro terms to a RefSeq GFF3, matched on NCBI GeneID.
    annotate-from-uniprot.py <uniprot.tsv> <in.gff> <out.gff>
A term goes on every gene and mRNA with that GeneID. Standard library only."""
import csv, re, sys
from collections import defaultdict
tsv, gff, out = sys.argv[1:4]
terms = defaultdict(set)
with open(tsv, newline="") as h:
    for r in csv.DictReader(h, delimiter="\t"):
        for g in filter(None, r["GeneID"].replace(" ", "").split(";")):
            terms[g].update(filter(None, r["InterPro"].replace(" ", "").split(";")))
            terms[g].update(filter(None, r["Gene Ontology IDs"].replace(" ", "").split(";")))
n = 0
with open(gff) as a, open(out, "w") as b:
    for line in a:
        f = line.rstrip("\n").split("\t")
        if len(f) == 9 and f[2] in ("gene", "mRNA"):
            m = re.search(r"GeneID:(\d+)", f[8])
            if m and terms.get(m.group(1)):
                f[8] += ";uniprot_terms=" + ",".join(sorted(terms[m.group(1)]))
                n += f[2] == "gene"
                line = "\t".join(f) + "\n"
        b.write(line)
print("genes annotated", n, "uniprot gene ids", len(terms))
