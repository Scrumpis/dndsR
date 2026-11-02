#!/usr/bin/env python3
"""Add public GO and InterPro terms to the gene lines of an A. thaliana RefSeq GFF3.

    annotate-gff.py <in.gff> <tair.gaf.gz> <uniprot.tsv> <out.gff> <truth.tsv>

GO comes from the GO Consortium's TAIR GAF, without NOT qualifiers or ND (no data).
InterPro comes from a UniProt proteome table with TAIR and InterPro columns.
truth.tsv records the expected terms per locus, to check append_annotations against.
Standard library only.
"""
import gzip
import re
import sys
from collections import defaultdict

gff_in, gaf, uniprot, gff_out, truth = sys.argv[1:6]
locus = re.compile(r"^AT[1-5CM]G\d{5}$")
go, ipr = defaultdict(set), defaultdict(set)

with gzip.open(gaf, "rt") as handle:
    for line in handle:
        if line.startswith("!"):
            continue
        f = line.rstrip("\n").split("\t")
        gene = f[1].upper()
        if locus.match(gene) and "NOT" not in f[3].split("|") and f[6] != "ND":
            go[gene].add(f[4])

with open(uniprot) as handle:
    head = handle.readline().rstrip("\n").split("\t")
    tair_col, ipr_col = head.index("TAIR"), head.index("InterPro")
    for line in handle:
        f = line.rstrip("\n").split("\t")
        terms = [x for x in f[ipr_col].split(";") if x]
        for gene in (x.upper() for x in f[tair_col].split(";") if x):
            if locus.match(gene):
                ipr[gene].update(terms)

seen = set()
with open(gff_in) as src, open(gff_out, "w") as out:
    for line in src:
        f = line.rstrip("\n").split("\t")
        if len(f) == 9 and f[2] == "gene":
            match = re.search(r"(?:^|;)locus_tag=([^;]+)", f[8])
            gene = match.group(1).upper() if match else ""
            extra = []
            if go.get(gene):
                extra.append("Ontology_term=" + ",".join(sorted(go[gene])))
            if ipr.get(gene):
                extra.append("interpro=" + ",".join(sorted(ipr[gene])))
            if extra:
                seen.add(gene)
                f[8] += ";" + ";".join(extra)
                line = "\t".join(f) + "\n"
        out.write(line)

with open(truth, "w") as out:
    out.write("gene\tgo\tipr\n")
    for gene in sorted(seen):
        out.write(f"{gene}\t{';'.join(sorted(go.get(gene, [])))}\t{';'.join(sorted(ipr.get(gene, [])))}\n")
print(f"GAF loci {len(go)}, UniProt loci {len(ipr)}, GFF genes annotated {len(seen)}")
