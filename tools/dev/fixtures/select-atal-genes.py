#!/usr/bin/env python3
"""Choose the genes of the A. thaliana - A. lyrata regression fixture.

Reads a full-genome dndsR run on the public RefSeq pair (outputs of extract_cds and
calculate_dnds) and the NCBI GFF3 files, and writes genes.tsv: one row per chosen gene,
with the stratum that chose it. Strata cover the situations the pipeline must handle:
dS from near zero to saturated, dN/dS above 1 and undefined, reliability flags, strand,
exon count, isoform count, short and long CDS, genes without a reciprocal partner, and
genes extract_cds excludes. Sampling is seeded, so the same inputs give the same list.
The list, not this script, defines the fixture: it is committed beside it.

  select-atal-genes.py <run dir with AtAl_*> <At gff.gz> <Al gff.gz> <out genes.tsv>

Standard library only; uses no dndsR code.
"""
import csv, gzip, random, sys
from collections import defaultdict

run, at_gff, al_gff, out = sys.argv[1:5]
random.seed(20261120)

def gff_genes(path):
    """gene id (ID without 'gene-') -> (seqid, start, end, strand, n_exons_max), for
    nuclear genes only: organelle sequences (region genome=chloroplast or mitochondrion)
    are left out, as docs/tutorials/public-genomes.md removes them before dndsR."""
    genes, tx_parent, exons, organelle = {}, {}, defaultdict(int), set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "region" and a.get("genome") in ("chloroplast", "mitochondrion"):
                organelle.add(f[0])
            if f[0] in organelle: continue
            if f[2] == "gene" and a.get("gene_biotype") == "protein_coding":
                genes[a["ID"][5:]] = [f[0], int(f[3]), int(f[4]), f[6], 0, a.get("Name", "")]
            elif f[2] == "mRNA" and a.get("Parent", "").startswith("gene-"):
                tx_parent[a["ID"]] = a["Parent"][5:]
            elif f[2] == "CDS" and "Parent" in a:
                exons[a["Parent"]] += 1
    for tx, g in tx_parent.items():
        if g in genes: genes[g][4] = max(genes[g][4], exons.get(tx, 0))
    return genes

def read(path):
    with open(path) as fh: return list(csv.DictReader(fh, delimiter="\t"))

at, al = gff_genes(at_gff), gff_genes(al_gff)
dnds = read(f"{run}/AtAl_dnds.tsv")
rel = {(r["query_id"], r["subject_id"]): r for r in read(f"{run}/AtAl_dnds_reliability.tsv")}
qc_at = {r["id"]: r for r in read(f"{run}/Athaliana_CDS_qc.tsv")}
qc_al = {r["id"]: r for r in read(f"{run}/Alyrata_CDS_qc.tsv")}

def num(x):
    try: return float(x)
    except ValueError: return float("nan")

pairs = [p for p in dnds if p["query_id"] in at and p["subject_id"] in al]
ds_values = sorted(num(p["dS"]) for p in pairs if num(p["dS"]) == num(p["dS"]))
median_ds = ds_values[len(ds_values) // 2]

strata = {
    "ds_below_0.02":   lambda p: num(p["dS"]) < 0.02,
    "ds_0.02_0.1":     lambda p: 0.02 <= num(p["dS"]) < 0.1,
    "ds_0.1_0.3":      lambda p: 0.1 <= num(p["dS"]) < 0.3,
    "ds_0.3_0.8":      lambda p: 0.3 <= num(p["dS"]) < 0.8,
    "ds_above_0.8":    lambda p: num(p["dS"]) >= 0.8,
    "ds_over_3x_median": lambda p: num(p["dS"]) > 3 * median_ds,
    "dnds_above_1":    lambda p: num(p["dNdS"]) > 1,
    "dnds_undefined":  lambda p: num(p["dNdS"]) != num(p["dNdS"]),
    "dn_zero":         lambda p: num(p["dN"]) == 0,
    "flag_few_codons": lambda p: "few_codons" in rel.get((p["query_id"], p["subject_id"]), {}).get("flags", ""),
    "flag_few_syn":    lambda p: "few_synonymous_differences" in rel.get((p["query_id"], p["subject_id"]), {}).get("flags", ""),
    "minus_strand":    lambda p: at[p["query_id"]][3] == "-",
    "single_exon":     lambda p: at[p["query_id"]][4] == 1,
    "many_exons":      lambda p: at[p["query_id"]][4] >= 10,
    "many_isoforms":   lambda p: int(qc_at.get(p["query_id"], {}).get("n_isoforms", "1") or 1) >= 4,
    "short_cds":       lambda p: int(qc_at.get(p["query_id"], {}).get("length_nt", "0") or 0) < 300,
    "long_cds":        lambda p: int(qc_at.get(p["query_id"], {}).get("length_nt", "0") or 0) > 4500,
}
per_stratum = {k: 25 for k in strata}
per_stratum.update({"ds_0.02_0.1": 40, "ds_0.1_0.3": 60, "ds_0.3_0.8": 40})

chosen = {}   # (species, gene) -> stratum
def take(species, gene, stratum):
    chosen.setdefault((species, gene), stratum)

for name, test in strata.items():
    pool = [p for p in pairs if test(p) and (("At", p["query_id"]) not in chosen)]
    for p in random.sample(pool, min(per_stratum[name], len(pool))):
        take("At", p["query_id"], name); take("Al", p["subject_id"], name)

# Gene symbols (Name) shared by two genes. Whole RefSeq genomes have them, and with
# them Bioconductor names genes by ID (locus tag), as in the full-genome run; a file in
# which every Name is unique is named by Name instead. Keeping a few makes the fixture
# take the whole-genome path. Both genes of a symbol are taken, with their partners.
by_name = defaultdict(list)
for p in pairs: by_name[at[p["query_id"]][5]].append(p)
dups = sorted(n for n, ps in by_name.items() if n and len(ps) == 2)
for n in random.sample(dups, min(3, len(dups))):
    for p in by_name[n]:
        chosen[("At", p["query_id"])] = "duplicate_gene_name"; take("Al", p["subject_id"], "duplicate_gene_name")

# Genes with a valid CDS but no reciprocal partner in the full-genome run.
paired_at = {p["query_id"] for p in dnds}; paired_al = {p["subject_id"] for p in dnds}
for sp, qc, genes, paired in (("At", qc_at, at, paired_at), ("Al", qc_al, al, paired_al)):
    pool = sorted(g for g, r in qc.items() if r.get("written") == "TRUE" and g not in paired and g in genes)
    for g in random.sample(pool, min(30, len(pool))): take(sp, g, "no_partner")
    # Genes extract_cds excludes: internal stop, or length not a multiple of 3.
    for col, label in (("internal_stop", "TRUE"), ("multiple_of_3", "FALSE")):
        pool = sorted(g for g, r in qc.items() if r.get(col) == label and g in genes)
        for g in random.sample(pool, min(12, len(pool))): take(sp, g, f"excluded_{col}")

with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["species", "gene_id", "stratum"])
    for (sp, g), s in sorted(chosen.items()): w.writerow([sp, g, s])
counts = defaultdict(int)
for (sp, _), s in chosen.items(): counts[sp] += 1
print(f"median dS {median_ds:.4f}; chosen: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
