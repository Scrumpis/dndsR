#!/usr/bin/env python3
"""Writes the stage 2B scenario and unit tables (docs/dev/phylo-calibration-stage2b.md).

  python3 -I phylo-stage2b-tables.py <out dir>

Every scenario has its own simulation seed (8001-8050), none shared with stage 2A
(7001-7135), so stage 2B's replicates are independent of those that produced 2A's findings.
Corrupted selection cells (B5) reuse the seed of their clean B3 cell, so the same genes are
compared with and without corruption. B2 is written twice, with 1 and with 3 starts; the
B1 audit rule (D4) chooses one before B2 is launched.
"""
import csv
import os
import sys

out = sys.argv[1]
S = []  # scenario, cell, n, seed, corrupt, corrupt_taxon, corrupt_seed, kind
def sc(name, cell, n, seed, kind, corrupt="none", taxon="-", cseed="-"):
    S.append([name, cell, n, seed, corrupt, taxon, cseed, kind])

# B1 audit (null-like), 300 genes. N0_L25 and N0_L50_div2 are also B2 grid cells.
sc("b_N0_L25", "N0_L25", 300, 8001, "null")
sc("b_N0_L50_div2", "N0_L50_div2", 300, 8002, "null")
sc("b_N1_L50", "N1_L50", 300, 8003, "null")
sc("b_M1_d02", "M1_mnp", 300, 8004, "misspec")
sc("b_R1_At", "N0_L50", 300, 8005, "corrupt", "internal", "At", 8105)
sc("b_M4", "M4_indel", 300, 8006, "misspec")
# B2 grid, the other 7 cells, 300 genes.
for i, c in enumerate(["N0_L25_div05", "N0_L25_div2", "N0_L50_div05", "N0_L50", "N0_L75_div05", "N0_L75", "N0_L75_div2"]):
    sc("b_" + c, c, 300, 8011 + i, "null")
# B3 power, 120 genes.
power = ["A2_p20", "A3_p05", "A3_p10", "A3_p20", "A5_p10", "A5_p20", "A10_p05", "A10_p10", "A10_p20"]
pseed = {c: 8021 + i for i, c in enumerate(power)}
for c in power:
    sc("b_" + c, c, 120, pseed[c], "power")
# B5 selection with corruption: the B3 genes, one internal segment in At.
for i, c in enumerate(["A3_p20", "A5_p10", "A10_p10", "A10_p20"]):
    sc("b_" + c + "_R1", c, 120, pseed[c], "corrupt_alt", "internal", "At", 8121 + i)
# B5b clustered selection, B4 multinucleotide, B6 indels with selection.
sc("b_AC10_p05", "AC10_p05", 120, 8031, "power")
sc("b_AC5_p10", "AC5_p10", 120, 8032, "power")
sc("b_M1_d03", "M1_mnp_d03", 150, 8041, "misspec")
sc("b_A10_p10_indel", "A10_p10_indel", 120, 8042, "power")
# B7 genome-like mixture: 720 null, 80 selected; per-taxon corruption at the P3A rates.
sc("b_MIX", "MIX:N0_L50:720+A3_p20:20+A5_p10:20+A10_p05:20+A10_p10:20", 800, 8050, "mixture",
   "pool", "At=0.032,Al=0.020,Cr=0.056,Es=0.083,Br=0.154", 8150)

U = []  # unit, block, scenario, config, filter, seeds, fit
def un(block, s, c, f="none", seeds="1,2,3", fit="all"):
    tag = {"1,2,3": "s123", "4,5,6": "s456", "1": "s1"}[seeds]
    U.append(["__".join([block, s, c, f, tag]), block, s, c, f, seeds, fit])

for s in ["b_N0_L25", "b_N0_L50_div2", "b_N1_L50", "b_M1_d02"]:
    for c in ["C1", "C3"]:
        un("B1", s, c)
for c in ["C1", "C3"]:
    un("B1", "b_R1_At", c, "mask")
for c in ["C1", "C3"]:
    un("B6", "b_M4", c)
un("B6", "b_M4", "C1", "true")
un("B6", "b_M4", "C1", "mask", fit="changed:none")
un("B6", "b_M4", "C1", "gapflank", fit="changed:none")
for c in power:
    s = "b_" + c
    for cf in ["C1", "C3"]:
        un("B3", s, cf)
        un("B3", s, cf, "mask", fit="changed:none")
    un("B3", s, "C1", "mask_fg", fit="changed:mask")
for c in ["A3_p20", "A10_p05", "A10_p10", "A10_p20"]:
    un("B3b", "b_" + c, "C0")
for c in ["A10_p10", "A10_p20"]:
    un("B4", "b_" + c, "C4")
    un("B4", "b_" + c, "C1", seeds="4,5,6")
for c in ["C1", "C4"]:
    un("B4", "b_M1_d03", c)
for c in ["A3_p20", "A5_p10", "A10_p10", "A10_p20"]:
    s = "b_" + c + "_R1"
    un("B5", s, "C1"); un("B5", s, "C1", "mask")
    un("B5", s, "C1", "mask_fg", fit="changed:mask")
    if c == "A10_p10":
        un("B5", s, "C3"); un("B5", s, "C3", "mask")
for s in ["b_AC10_p05", "b_AC5_p10"]:
    un("B5b", s, "C1"); un("B5b", s, "C1", "mask", fit="changed:none")
    un("B5b", s, "C1", "mask_fg", fit="changed:mask"); un("B5b", s, "C3")
un("B6", "b_A10_p10_indel", "C1")
un("B6", "b_A10_p10_indel", "C1", "mask", fit="changed:none")
un("B6", "b_A10_p10_indel", "C1", "gapflank", fit="changed:none")
for c in ["C1", "C3"]:
    un("B7", "b_MIX", c)
    un("B7", "b_MIX", c, "mask", fit="changed:none")
    un("B7", "b_MIX", c, "mask_fg", fit="changed:mask")

# Per-unit cost estimates (core-hours), from stage 2A's measured core-seconds per fit, for the
# budget guard in phylo-stage2-run.sb. Filtered units fit only changed alignments; the shares
# changed are assumptions from 2A's flag rates and are replaced by measured costs as units finish.
COST = {("C0", "power"): 97, ("C1", "null"): 78, ("C1", "power"): 169, ("C1", "misspec"): 111, ("C1", "corrupt_alt"): 170,
        ("C3", "null"): 111, ("C3", "power"): 235, ("C3", "misspec"): 167, ("C3", "corrupt_alt"): 280,
        ("C4", "power"): 318, ("C4", "misspec"): 135}
LENF = {"N0_L25": 0.67, "N0_L25_div05": 0.67, "N0_L25_div2": 0.67, "N0_L75": 1.43, "N0_L75_div05": 1.43, "N0_L75_div2": 1.43}
SC = {r[0]: r for r in S}
def estimate(u):
    _, _, s, c, f, seeds, fit = u
    name, cell, n, _, _, _, _, kind = SC[s]
    if kind == "mixture":
        per = {"C1": 0.9 * 78 + 0.1 * 169, "C3": 0.9 * 111 + 0.1 * 235}[c]
    else:
        k = "null" if kind == "corrupt" else kind
        per = COST.get((c, k)) or COST[(c, "misspec")]
    share = 1.0
    if fit != "all":
        if s == "b_MIX": share = 0.30 if f == "mask" else 0.25
        elif s.startswith("b_AC"): share = 0.5 if f == "mask" else 0.05
        elif f == "gapflank": share = 0.95
        elif s in ("b_M4", "b_A10_p10_indel"): share = 0.08
        else: share = 0.03 if f == "mask" else 0.02
    return round(n * len(seeds.split(",")) * per * LENF.get(cell, 1) * share / 3600, 2)

EST = []
def write(name, rows, header):
    with open(os.path.join(out, name), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header); w.writerows(rows)
    if name.startswith("phylo-stage2b-units"):
        EST.extend([r[0], estimate(r)] for r in rows)

write("phylo-stage2b-scenarios.tsv", S, ["scenario", "cell", "n", "seed", "corrupt", "corrupt_taxon", "corrupt_seed", "kind"])
hdr = ["unit", "block", "scenario", "config", "filter", "seeds", "fit"]
write("phylo-stage2b-units.tsv", U, hdr)
grid = ["b_N0_L25_div05", "b_N0_L25_div2", "b_N0_L50_div05", "b_N0_L50", "b_N0_L75_div05", "b_N0_L75", "b_N0_L75_div2"]
for seeds, tag in [("1", "1start"), ("1,2,3", "3start")]:
    U = []
    for s in grid:
        for c in ["C1", "C3"]:
            un("B2", s, c, seeds=seeds)
    write("phylo-stage2b-units-B2-%s.tsv" % tag, U, hdr)
write("phylo-stage2b-estimates.tsv", EST, ["unit", "core_hours"])
