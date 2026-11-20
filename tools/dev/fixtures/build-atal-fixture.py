#!/usr/bin/env python3
"""Build the A. thaliana - A. lyrata regression fixture from NCBI RefSeq files.

For each gene in genes.tsv, takes its genomic span plus FLANK bases on each side,
merges windows that overlap, and writes each merged window as one contig of real
sequence. Every GFF3 feature whose top-level gene lies wholly inside a window is kept,
with coordinates shifted to the contig; genes cut by a window edge are dropped whole.
So neighbouring genes, isoforms, UTRs, non-coding genes, and pseudogenes inside a
window stay as annotated. Contigs are named <seqid>_<start>_<end>, 1-based inclusive
coordinates on the source sequence. Output is gzipped, written with fixed timestamps
so a rebuild is byte-identical.

  build-atal-fixture.py <genes.tsv> <At fna.gz> <At gff.gz> <Al fna.gz> <Al gff.gz> <out dir>

Standard library only; uses no dndsR code.
"""
import csv, gzip, io, os, sys
from collections import defaultdict

FLANK = 1000
genes_tsv, at_fna, at_gff, al_fna, al_gff, out = sys.argv[1:7]
os.makedirs(out, exist_ok=True)

def read_gff(path):
    lines = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f = line.rstrip("\n").split("\t")
            if len(f) == 9: lines.append(f)
    return lines

def attrs(f): return dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)

def root_of(lines):
    """feature index -> ID of its top-level ancestor (itself when it has no Parent)."""
    parent = {}
    ids = []
    for f in lines:
        a = attrs(f); ids.append(a.get("ID")); 
        if a.get("ID") is not None: parent[a["ID"]] = a.get("Parent", "").split(",")[0] or None
    def root(i):
        seen = 0
        while parent.get(i) and seen < 20: i = parent[i]; seen += 1
        return i
    roots = []
    for f, i in zip(lines, ids):
        if i is None: i = attrs(f).get("Parent", "").split(",")[0] or None
        roots.append(root(i) if i else None)
    return roots

def read_fasta(path, wanted):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name in wanted: seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name in wanted: buf.append(line.strip())
    if name in wanted: seqs[name] = "".join(buf)
    return seqs

def gz_write(path, text):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            gz.write(text.encode())

wanted = defaultdict(set)
with open(genes_tsv) as fh:
    for r in csv.DictReader(fh, delimiter="\t"): wanted[r["species"]].add(r["gene_id"])

summary = []
for sp, fna, gff in (("At", at_fna, at_gff), ("Al", al_fna, al_gff)):
    lines = read_gff(gff); roots = root_of(lines)
    # A gene can span several lines with one ID (multi-part genes): take their union.
    # Organelle sequences are never windowed (see select-atal-genes.py).
    organelle = {f[0] for f in lines if f[2] == "region" and attrs(f).get("genome") in ("chloroplast", "mitochondrion")}
    span = {}
    for f in lines:
        a = attrs(f)
        if f[2] in ("gene", "pseudogene") and a.get("ID") and f[0] not in organelle:
            s0, a0, b0 = span.get(a["ID"], (f[0], int(f[3]), int(f[4])))
            span[a["ID"]] = (f[0], min(a0, int(f[3])), max(b0, int(f[4]))) if s0 == f[0] else (s0, -1, -1)
    missing = [g for g in wanted[sp] if "gene-" + g not in span]
    if missing: sys.exit(f"{sp}: genes not in GFF: {missing[:5]}")
    seq_ids = {span["gene-" + g][0] for g in wanted[sp]}
    seqs = read_fasta(fna, seq_ids)
    by_seq = defaultdict(list)
    for g in wanted[sp]:
        s, a, b = span["gene-" + g]
        by_seq[s].append([max(1, a - FLANK), min(len(seqs[s]), b + FLANK)])
    windows = []
    for s, iv in by_seq.items():
        iv.sort()
        cur = iv[0]
        for a, b in iv[1:]:
            if a <= cur[1] + 1: cur[1] = max(cur[1], b)
            else: windows.append((s, *cur)); cur = [a, b]
        windows.append((s, *cur))
    windows.sort(key=lambda w: (w[0], w[1]))
    win_of = defaultdict(list)
    for w in windows: win_of[w[0]].append(w)
    def window_for(s, a, b):
        for w in win_of.get(s, ()):
            if w[1] <= a and b <= w[2]: return w
        return None
    keep_root = {}
    for rid, (s, a, b) in span.items():
        w = window_for(s, a, b)
        if w: keep_root[rid] = w
    fa, gf = io.StringIO(), io.StringIO()
    gf.write("##gff-version 3\n")
    for s, a, b in windows:
        name = f"{s}_{a}_{b}"
        gf.write(f"##sequence-region {name} 1 {b - a + 1}\n")
        sub = seqs[s][a - 1:b]
        fa.write(f">{name}\n" + "\n".join(sub[i:i + 80] for i in range(0, len(sub), 80)) + "\n")
    n_feat, kept_genes = 0, set()
    for f, r in zip(lines, roots):
        w = keep_root.get(r)
        if not w: continue
        s, a, b = w
        if not (a <= int(f[3]) and int(f[4]) <= b):
            sys.exit(f"{sp}: feature outside its window: {f[0]} {f[2]} {f[3]}-{f[4]}")
        g = list(f); g[0] = f"{s}_{a}_{b}"; g[3] = str(int(f[3]) - a + 1); g[4] = str(int(f[4]) - a + 1)
        gf.write("\t".join(g) + "\n"); n_feat += 1
        if f[2] == "gene": kept_genes.add(r)
    sp_name = {"At": "Athaliana", "Al": "Alyrata"}[sp]
    gz_write(os.path.join(out, f"{sp_name}.fasta.gz"), fa.getvalue())
    gz_write(os.path.join(out, f"{sp_name}.gff3.gz"), gf.getvalue())
    bases = sum(b - a + 1 for _, a, b in windows)
    summary.append(f"{sp_name}: {len(wanted[sp])} chosen genes, {len(kept_genes)} genes kept in "
                   f"{len(windows)} contigs, {bases} bases, {n_feat} features")
print("\n".join(summary))
