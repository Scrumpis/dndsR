#!/usr/bin/env python3
"""Counts transcripts whose first CDS segment has a nonzero GFF3 phase, which marks a partial gene.
dndsR does not trim these to frame, so they can give wrong codons if they pass its stop check.

    check-first-exon-phase.py <annotation.gff3>

Standard library only.
"""
import sys

first = {}  # transcript: (position on the transcript, phase)
with open(sys.argv[1]) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9 or f[2] != "CDS":
            continue
        attrs = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        order = int(f[3]) if f[6] == "+" else -int(f[4])
        for parent in attrs.get("Parent", attrs.get("ID", "")).split(","):
            if parent not in first or order < first[parent][0]:
                first[parent] = (order, f[7])
partial = [t for t, (_, phase) in first.items() if phase not in ("0", ".")]
print(f"{len(partial)} of {len(first)} transcripts have a first CDS segment of nonzero phase")
for t in sorted(partial)[:20]:
    print(" ", t)
