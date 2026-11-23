#!/usr/bin/env python3
"""Add GO and InterPro terms from a UniProt reference proteome to the gene lines of a RefSeq GFF3.

    annotate-gff-uniprot.py <in.gff> <protein.faa> <uniprot.tsv> <out.gff> <report.tsv>

uniprot.tsv: UniProt REST columns Entry, RefSeq, GeneID, Gene Ontology IDs, InterPro, Sequence.
A UniProt entry is assigned to RefSeq genes, in this order, by:
  1. NCBI GeneID cross-reference (matched to the gene line's Dbxref GeneID);
  2. RefSeq protein cross-reference (matched through the CDS lines' protein_id to their gene);
  3. exact protein sequence (matched to protein.faa, then through protein_id to the gene).
Terms are written as `Ontology_term=` and `interpro=` (append_annotations reads both), as
annotate-gff.py does for A. thaliana. Every term on an entry goes to every gene it maps to.
report.tsv records how many entries and genes each route mapped, for the study's coverage gate.
Standard library only.
"""
import re
import sys
from collections import defaultdict

gff_in, faa, uniprot, gff_out, report = sys.argv[1:6]


def attrs(field):
    out = {}
    for part in field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


# Genes, mRNA -> gene, protein -> gene, GeneID -> gene
gene_lines, rna_gene, prot_gene, geneid_gene = {}, {}, {}, {}
with open(gff_in) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            continue
        a = attrs(f[8])
        if f[2] == "gene" and a.get("gene_biotype") == "protein_coding":
            gid = a.get("ID")
            gene_lines[gid] = True
            m = re.search(r"GeneID:(\d+)", a.get("Dbxref", ""))
            if m:
                geneid_gene[m.group(1)] = gid
        elif f[2] in ("mRNA", "transcript") and "Parent" in a:
            rna_gene[a["ID"]] = a["Parent"]
        elif f[2] == "CDS" and "Parent" in a:
            pid = a.get("protein_id") or a.get("Name") or re.sub(r"^cds-", "", a.get("ID", ""))
            if pid:
                prot_gene[pid] = rna_gene.get(a["Parent"], a["Parent"])

seq_prot = defaultdict(list)
with open(faa) as fh:
    pid, buf = None, []
    for line in fh:
        if line.startswith(">"):
            if pid:
                seq_prot["".join(buf)].append(pid)
            pid, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if pid:
        seq_prot["".join(buf)].append(pid)

go, ipr = defaultdict(set), defaultdict(set)
counts = defaultdict(int)
with open(uniprot) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    col = {h: i for i, h in enumerate(head)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        counts["entries"] += 1
        terms_go = [t.strip() for t in f[col["Gene Ontology IDs"]].split(";") if t.strip()]
        terms_ipr = [t.strip() for t in f[col["InterPro"]].split(";") if t.strip()]
        if not terms_go and not terms_ipr:
            counts["entries_without_terms"] += 1
            continue
        genes, route = set(), None
        for g in (x.strip() for x in f[col["GeneID"]].split(";") if x.strip()):
            if g in geneid_gene:
                genes.add(geneid_gene[g])
        if genes:
            route = "geneid"
        else:
            for p in (x.strip().split(" ")[0] for x in f[col["RefSeq"]].split(";") if x.strip()):
                if p in prot_gene:
                    genes.add(prot_gene[p])
            if genes:
                route = "refseq_protein"
            else:
                for p in seq_prot.get(f[col["Sequence"]].strip().upper(), []):
                    if p in prot_gene:
                        genes.add(prot_gene[p])
                if genes:
                    route = "sequence"
        if not genes:
            counts["entries_unmapped"] += 1
            continue
        counts["entries_by_" + route] += 1
        for g in genes:
            go[g].update(terms_go)
            ipr[g].update(terms_ipr)

annotated = 0
with open(gff_in) as src, open(gff_out, "w") as out:
    for line in src:
        f = line.rstrip("\n").split("\t")
        if len(f) == 9 and f[2] == "gene":
            gid = attrs(f[8]).get("ID")
            extra = []
            if go.get(gid):
                extra.append("Ontology_term=" + ",".join(sorted(go[gid])))
            if ipr.get(gid):
                extra.append("interpro=" + ",".join(sorted(ipr[gid])))
            if extra:
                annotated += 1
                f[8] += ";" + ";".join(extra)
                line = "\t".join(f) + "\n"
        out.write(line)

counts["protein_coding_genes"] = len(gene_lines)
counts["genes_with_go"] = sum(1 for g in gene_lines if go.get(g))
counts["genes_with_ipr"] = sum(1 for g in gene_lines if ipr.get(g))
counts["genes_annotated"] = annotated
with open(report, "w") as out:
    out.write("measure\tvalue\n")
    for k in sorted(counts):
        out.write(f"{k}\t{counts[k]}\n")
print("\n".join(f"{k}\t{counts[k]}" for k in sorted(counts)))
