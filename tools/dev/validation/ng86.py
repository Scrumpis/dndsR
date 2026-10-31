#!/usr/bin/env python3
"""Independent Nei-Gojobori (1986) dN/dS for pairs of aligned CDS. Standard library only.

Shares no code with dndsR or orthologr. Used to establish expected values.
  ng86.py query.fasta subject.fasta   (pairs are matched by the name after the first '_')
Terminal stop codons are excluded. Pairs with internal stops, unequal length, or a length
not divisible by 3 are reported as such rather than computed.
"""
import itertools
import math
import sys

BASES = "TCAG"
CODE = dict(zip(("".join(c) for c in itertools.product(BASES, repeat=3)),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def site_counts(codon):
    """Synonymous sites in a codon; changes to a stop codon are excluded."""
    s = 0.0
    for pos in range(3):
        syn = total = 0
        for base in BASES:
            if base == codon[pos]:
                continue
            mutant = codon[:pos] + base + codon[pos + 1:]
            if CODE[mutant] == "*":
                continue
            total += 1
            syn += CODE[mutant] == CODE[codon]
        if total:
            s += syn / total
    return s


def differences(c1, c2):
    """Mean synonymous and nonsynonymous differences over all shortest paths avoiding stops."""
    diff = [i for i in range(3) if c1[i] != c2[i]]
    if not diff:
        return 0.0, 0.0
    sd = nd = paths = 0
    for order in itertools.permutations(diff):
        cur, s, n, ok = c1, 0, 0, True
        for pos in order:
            nxt = cur[:pos] + c2[pos] + cur[pos + 1:]
            if CODE[nxt] == "*":
                ok = False
                break
            if CODE[nxt] == CODE[cur]:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            sd, nd, paths = sd + s, nd + n, paths + 1
    return (sd / paths, nd / paths) if paths else (0.0, 0.0)


def jukes_cantor(p):
    return math.inf if p >= 0.75 else -0.75 * math.log(1 - 4 * p / 3)


def ng86(a, b):
    if len(a) != len(b):
        return "unequal_length"
    if len(a) % 3:
        return "length_not_multiple_of_3"
    cod_a = [a[i:i + 3] for i in range(0, len(a), 3)]
    cod_b = [b[i:i + 3] for i in range(0, len(b), 3)]
    if CODE[cod_a[-1]] == "*" and CODE[cod_b[-1]] == "*":
        cod_a, cod_b = cod_a[:-1], cod_b[:-1]
    if any(CODE[c] == "*" for c in cod_a + cod_b):
        return "internal_stop"
    sites_s = sum((site_counts(x) + site_counts(y)) / 2 for x, y in zip(cod_a, cod_b))
    sites_n = 3 * len(cod_a) - sites_s
    sd = nd = 0.0
    for x, y in zip(cod_a, cod_b):
        s, n = differences(x, y)
        sd, nd = sd + s, nd + n
    dn, ds = jukes_cantor(nd / sites_n), jukes_cantor(sd / sites_s)
    return sd, nd, sites_s, sites_n, dn, ds


def main(query, subject):
    q, s = read_fasta(query), read_fasta(subject)
    key = lambda n: n.split("_", 1)[1] if "_" in n else n
    subjects = {key(n): v for n, v in s.items()}
    print("pair\tSd\tNd\tS\tN\tdN\tdS\tdNdS\tstatus")
    for name, seq in q.items():
        if key(name) not in subjects:
            continue
        res = ng86(seq, subjects[key(name)])
        if isinstance(res, str):
            print(f"{key(name)}\t\t\t\t\t\t\t\t{res}")
            continue
        sd, nd, ss, nn, dn, ds = res
        ratio = f"{dn / ds:.6f}" if ds > 0 else "NA"
        print(f"{key(name)}\t{sd:.1f}\t{nd:.1f}\t{ss:.4f}\t{nn:.4f}\t{dn:.6f}\t{ds:.6f}\t{ratio}\tok")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
