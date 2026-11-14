#!/usr/bin/env python3
"""Locally dissimilar windows in codon alignments: profile, filter, and corrupt.

Independent of dndsR's code (base Python). The window rule is the one documented for
dndsR's `dissimilar_window` flag (R/phylo_alignment.R, .phylo_qc_limits): 10-codon
windows with at least 8 determined codons, of which at least 70% have an amino acid
unlike every other sequence at that column. Gaps, ambiguous codons, and a final stop
count as undetermined.

  phylo-windows.py profile <alignment dir> <out.tsv>
      One row per sequence: windows (count, comparable with dndsR's dissimilar_windows),
      and its merged segments (start, end, 1-based codons) with whether each touches an
      end of the sequence and whether a gap in that sequence lies within 3 codons.
  phylo-windows.py filter <in dir> <out dir> mask|mask+gappy
      mask: codons of a sequence inside its dissimilar windows become gaps.
      mask+gappy: also drop columns gapped in more than half the sequences.
      mask_fg <taxa>: as mask, in the sequences of the comma-separated taxa only (stage 2B).
      gapflank: as mask, and also, in each sequence, the codons within 3 of a gap run inside
      that sequence (not at its ends), where an indel may have been misplaced (stage 2B).
      Every mode reads the alignment only. Writes <out dir>/filter_log.tsv (codons masked,
      columns dropped, columns kept, codons present before filtering, per alignment).
  phylo-windows.py corrupt-pool <in dir> <out dir> <rates> <segments.tsv> <seed>
      Stage 2B genome-like corruption: each sequence independently, with its taxon's
      probability from <rates> (e.g. At=0.032,Al=0.02), gets one segment: terminal (18%),
      frameshift (48%), or internal (34%), the shares of the P3A foreground segments that
      touch an end, are internal next to a gap, and are internal otherwise.
      Writes <out dir>/corrupt_log.tsv.
  phylo-windows.py corrupt <in dir> <out dir> <kind> <taxon> <segments.tsv> <seed>
      kind: internal, terminal, or frameshift. Replaces one segment of <taxon> per
      alignment. Segment lengths are drawn from <segments.tsv> (a profile's foreground
      segments). internal and terminal use the same taxon's codons from another
      alignment in the folder (realistic composition, unrelated sequence); frameshift
      reads the segment one nucleotide out of frame, as compensating indels would.
      Writes <out dir>/corrupt_log.tsv.
"""
import csv
import os
import random
import sys

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
W, W_MIN, W_MISMATCH = 10, 8, 0.7


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            names.append(line[1:].split()[0]); seqs.append([])
        else:
            seqs[-1].append(line.upper())
    return names, ["".join(s) for s in seqs]


def write_fasta(path, names, seqs):
    with open(path, "w") as f:
        for n, s in zip(names, seqs):
            f.write(">" + n + "\n" + s + "\n")


def codons(s):
    return [s[i:i + 3] for i in range(0, len(s), 3)]


def translate(cods):
    aa = [CODE.get(c) for c in cods]
    # A stop ending the sequence is undetermined (dndsR masks it); internal stops do not occur.
    last = max((k for k, c in enumerate(cods) if c != "---"), default=None)
    if last is not None and aa[last] == "*":
        aa[last] = None
    return aa


def windows(seqs):
    """Per sequence: list of window starts (0-based) that are dissimilar."""
    aa = [translate(codons(s)) for s in seqs]
    n, L = len(aa), len(aa[0])
    out = []
    for i in range(n):
        unlike = []
        for k in range(L):
            a = aa[i][k]
            others = [aa[j][k] for j in range(n) if j != i and aa[j][k] is not None]
            unlike.append(a is not None and len(others) > 0 and a not in others)
        starts = []
        for s in range(0, L - W + 1):
            k = sum(1 for t in range(s, s + W) if aa[i][t] is not None)
            if k >= W_MIN and sum(unlike[s:s + W]) / k >= W_MISMATCH:
                starts.append(s)
        out.append(starts)
    return out


def segments(starts):
    segs = []
    for s in starts:
        if segs and s <= segs[-1][1] + 1:
            segs[-1][1] = max(segs[-1][1], s + W - 1)
        else:
            segs.append([s, s + W - 1])
    return segs


def taxon(name):
    return name.split("|")[0]


def profile(adir, out):
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["alignment", "taxon", "seq_name", "codons", "windows", "segment_start", "segment_end",
                    "touches_end", "gap_within_3"])
        for fn in sorted(os.listdir(adir)):
            if not fn.endswith((".fasta", ".fas")):
                continue
            names, seqs = read_fasta(os.path.join(adir, fn))
            ws = windows(seqs)
            for n, s, st in zip(names, seqs, ws):
                cods = codons(s); L = len(cods)
                det = [k for k, c in enumerate(cods) if c != "---"]
                first, last = (det[0], det[-1]) if det else (0, L - 1)
                segs = segments(st)
                if not segs:
                    w.writerow([fn, taxon(n), n, L, 0, "", "", "", ""]); continue
                for a, b in segs:
                    touches = a <= first + 2 or b >= last - 2
                    gap = any(cods[k] == "---" for k in range(max(0, a - 3), min(L, b + 4)) if first <= k <= last)
                    w.writerow([fn, taxon(n), n, L, len(st), a + 1, b + 1, int(touches), int(gap)])


def gap_flanks(cods, d=3):
    """Codon positions within d of a gap run inside the sequence (not at its ends)."""
    det = [k for k, c in enumerate(cods) if c != "---"]
    if not det:
        return set()
    first, last, L, out, k = det[0], det[-1], len(cods), set(), det[0]
    while k <= last:
        if cods[k] == "---":
            e = k
            while e + 1 <= last and cods[e + 1] == "---":
                e += 1
            out.update(j for j in range(max(0, k - d), k))
            out.update(j for j in range(e + 1, min(L, e + 1 + d)))
            k = e + 1
        else:
            k += 1
    return out


def filt(indir, outdir, mode, fg=None):
    os.makedirs(outdir, exist_ok=True)
    fg = set(fg.split(",")) if fg else None
    log = [["alignment", "codons_masked", "columns_dropped", "columns_kept", "codons_present"]]
    for fn in sorted(os.listdir(indir)):
        if not fn.endswith((".fasta", ".fas")) or fn.endswith(".true.fasta"):
            continue
        names, seqs = read_fasta(os.path.join(indir, fn))
        cods = [codons(s) for s in seqs]
        present = sum(c != "---" for cs in cods for c in cs)
        masked = 0
        flank = [gap_flanks(c) for c in cods] if mode == "gapflank" else None
        for i, st in enumerate(windows(seqs)):
            if mode == "mask_fg" and taxon(names[i]) not in fg:
                continue
            for a, b in segments(st):
                for k in range(a, b + 1):
                    if cods[i][k] != "---":
                        cods[i][k] = "---"; masked += 1
        if flank:
            for i, fl in enumerate(flank):
                for k in fl:
                    if cods[i][k] != "---":
                        cods[i][k] = "---"; masked += 1
        L = len(cods[0]); keep = list(range(L))
        if mode == "mask+gappy":
            keep = [k for k in range(L) if sum(c[k] == "---" for c in cods) <= len(cods) / 2]
        log.append([fn, masked, L - len(keep), len(keep), present])
        write_fasta(os.path.join(outdir, fn), names, ["".join(c[k] for k in keep) for c in cods])
    with open(os.path.join(outdir, "filter_log.tsv"), "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(log)


def corrupt(indir, outdir, kind, tx, seg_file, seed):
    rng = random.Random(int(seed))
    lengths = [int(r["segment_end"]) - int(r["segment_start"]) + 1 for r in csv.DictReader(open(seg_file), delimiter="\t")
               if r["segment_start"]]
    os.makedirs(outdir, exist_ok=True)
    files = sorted(fn for fn in os.listdir(indir) if fn.endswith(".fasta") and not fn.endswith(".true.fasta"))
    data = {fn: read_fasta(os.path.join(indir, fn)) for fn in files}
    log = [["alignment", "kind", "taxon", "start", "end", "donor"]]
    for idx, fn in enumerate(files):
        names, seqs = data[fn]
        i = [taxon(n) for n in names].index(tx)
        cods = codons(seqs[i]); L = len(cods)
        m = min(rng.choice(lengths), L // 2)
        if kind == "terminal":
            a = 0 if rng.random() < 0.5 else L - m
        else:
            a = rng.randrange(0, L - m + 1)
        b = a + m  # exclusive
        donor = ""
        if kind in ("internal", "terminal"):
            donor = files[(idx + 1) % len(files)]
            dn, ds = data[donor]
            dc = codons(ds[[taxon(n) for n in dn].index(tx)])
            off = rng.randrange(0, len(dc) - m + 1)
            new = dc[off:off + m]
        elif kind == "frameshift":
            nt = "".join(cods[a:b])
            shifted = nt[1:] + rng.choice("ACGT")
            new = codons(shifted)
            # A stop read out of frame would end the gene; replace it by a sense codon, as an
            # annotation that survived would have.
            new = [c if CODE[c] != "*" else rng.choice([k for k, v in CODE.items() if v != "*"]) for c in new]
        else:
            raise SystemExit("unknown kind " + kind)
        cods[a:b] = new
        seqs = list(seqs); seqs[i] = "".join(cods)
        write_fasta(os.path.join(outdir, fn), names, seqs)
        log.append([fn, kind, tx, a + 1, b, donor])
    with open(os.path.join(outdir, "corrupt_log.tsv"), "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(log)


def corrupt_pool(indir, outdir, rates, seg_file, seed):
    rng = random.Random(int(seed))
    rate = {k: float(v) for k, v in (x.split("=") for x in rates.split(","))}
    lengths = [int(r["segment_end"]) - int(r["segment_start"]) + 1 for r in csv.DictReader(open(seg_file), delimiter="\t")
               if r["segment_start"]]
    os.makedirs(outdir, exist_ok=True)
    files = sorted(fn for fn in os.listdir(indir) if fn.endswith(".fasta") and not fn.endswith(".true.fasta"))
    data = {fn: read_fasta(os.path.join(indir, fn)) for fn in files}
    sense = [k for k, v in CODE.items() if v != "*"]
    log = [["alignment", "kind", "taxon", "start", "end", "donor"]]
    for fn in files:
        names, seqs = data[fn]
        seqs = list(seqs)
        for i, n in enumerate(names):
            tx = taxon(n)
            if rng.random() >= rate.get(tx, 0.0):
                continue
            u = rng.random()
            kind = "terminal" if u < 0.18 else ("frameshift" if u < 0.66 else "internal")
            cods = codons(seqs[i]); L = len(cods)
            m = min(rng.choice(lengths), L // 2)
            a = (0 if rng.random() < 0.5 else L - m) if kind == "terminal" else rng.randrange(0, L - m + 1)
            b = a + m
            donor = ""
            if kind == "frameshift":
                nt = "".join(cods[a:b])
                new = codons(nt[1:] + rng.choice("ACGT"))
                new = [c if CODE.get(c, "X") != "*" else rng.choice(sense) for c in new]
            else:
                donor = rng.choice([f for f in files if f != fn])
                dn, ds = data[donor]
                dc = codons(ds[[taxon(x) for x in dn].index(tx)])
                off = rng.randrange(0, len(dc) - m + 1)
                new = dc[off:off + m]
            # Keep gaps where the sequence had them, as a real misassembled segment would not
            # fill an alignment gap; replace determined codons only.
            new = [c if old != "---" else "---" for c, old in zip(new, cods[a:b])]
            cods[a:b] = new
            seqs[i] = "".join(cods)
            log.append([fn, kind, tx, a + 1, b, donor])
        write_fasta(os.path.join(outdir, fn), names, seqs)
    with open(os.path.join(outdir, "corrupt_log.tsv"), "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(log)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "profile" and len(a) == 3:
        profile(a[1], a[2])
    elif a and a[0] == "filter" and len(a) == 4 and a[3] in ("mask", "mask+gappy", "gapflank"):
        filt(a[1], a[2], a[3])
    elif a and a[0] == "filter" and len(a) == 5 and a[3] == "mask_fg":
        filt(a[1], a[2], a[3], a[4])
    elif a and a[0] == "corrupt-pool" and len(a) == 6:
        corrupt_pool(*a[1:])
    elif a and a[0] == "corrupt" and len(a) == 7:
        corrupt(*a[1:])
    else:
        sys.exit(__doc__)
