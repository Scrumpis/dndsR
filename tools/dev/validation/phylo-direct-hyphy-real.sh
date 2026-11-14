#!/usr/bin/env bash
# Runs HyPhy BUSTED directly, without dndsR, on chosen orthogroups of a real phylo_busted
# run, and compares with dndsR's results. Trees are written here by hand for each set of taxa
# (species tree (((At,Al),Cr),(Es,Br)); foreground At, Al, and their stem), not taken from
# dndsR. Alignments are phylo_prepare's, with sequence names cut to the taxon by awk.
# Run on the same processor class as the dndsR run (HyPhy differs slightly between AVX-512
# and AVX2), inside the same image.
#   tools/dev/validation/phylo-direct-hyphy-real.sh <hyphy> <prepared dir or prepared.tar.gz> \
#     <busted_results.tsv> <orthogroups.txt> <work dir> [jobs]
# Writes <work dir>/direct_vs_dndsr.tsv.
set -euo pipefail
hyphy=$1; prepared=$2; results=$3; ogs=$4; work=$5; jobs=${6:-8}
mkdir -p "$work/aln"
if [[ $prepared == *.tar.gz ]]; then
  tar -xzf "$prepared" -C "$work/aln" --wildcards $(sed 's|.*|*/alignments/&.fasta|' "$ogs")
  adir=$(dirname "$(find "$work/aln" -name '*.fasta' | head -1)")
else
  adir=$prepared/alignments
fi
version=$("$hyphy" --version | grep -o '[0-9]*\.[0-9]*\.[0-9]*' | head -1)

tree_for() {   # hand-written, unrooted, for each set of taxa that occurs; anything else is an error
  case "$(printf '%s\n' "$@" | sort | tr '\n' ' ')" in
    "Al At Br Cr Es ") echo "((At{FG},Al{FG}){FG},Cr,(Es,Br));" ;;
    "At Br Cr Es ")    echo "(At{FG},Cr,(Es,Br));" ;;
    "Al Br Cr Es ")    echo "(Al{FG},Cr,(Es,Br));" ;;
    "Al At Cr Es ")    echo "((At{FG},Al{FG}){FG},Cr,Es);" ;;
    "Al At Br Cr ")    echo "((At{FG},Al{FG}){FG},Cr,Br);" ;;
    *) echo "unexpected taxa: $*" >&2; return 1 ;;
  esac
}

run_one() {
  og=$1; d="$work/$og"; mkdir -p "$d"
  awk '/^>/ { split(substr($0, 2), a, "|"); print ">" a[1]; next } { print }' "$adir/$og.fasta" > "$d/alignment.fas"
  tree_for $(grep '^>' "$d/alignment.fas" | cut -c2-) > "$d/tree.nwk"
  ( cd "$d" && "$hyphy" CPU=1 ENV="RANDOM_SEED=1" busted --alignment alignment.fas --tree tree.nwk \
      --branches FG --code Universal --srv Yes --rates 3 --syn-rates 3 --multiple-hits None --error-sink No \
      --grid-size 250 --starting-points 1 --kill-zero-lengths Yes --output busted.json > hyphy.log 2>&1 ) || true
}
export -f run_one tree_for; export hyphy work adir
xargs -P "$jobs" -I{} bash -c 'run_one {}' < "$ogs"

python3 -I - "$work" "$results" "$ogs" "$version" <<'EOF'
import csv, json, os, sys
work, results, ogs, version = sys.argv[1:]
want = [l.strip() for l in open(ogs) if l.strip()]
dn = {r["orthogroup"]: r for r in csv.DictReader(open(results), delimiter="\t") if r["orthogroup"] in want}
rows = ["orthogroup\thyphy_version\tdirect_status\tdirect_lrt\tdndsr_status\tdndsr_lrt\tabs_diff\tdirect_lnl_unconstrained\tdndsr_lnl_unconstrained\tdirect_tree"]
for og in want:
    j = os.path.join(work, og, "busted.json"); tree = open(os.path.join(work, og, "tree.nwk")).read().strip()
    try:   # HyPhy leaves an empty or partial JSON when it stops with an error
        x = json.load(open(j)); lrt = x["test results"]["LRT"]; lnl = x["fits"]["Unconstrained model"]["Log Likelihood"]; st = "tested"
    except (OSError, ValueError, KeyError):
        lrt = lnl = None; st = "hyphy_error"
    r = dn.get(og, {}); dl = r.get("lrt", "NA")
    diff = abs(float(dl) - lrt) if lrt is not None and dl not in ("", "NA") else None
    rows.append("\t".join(map(str, [og, version, st, "NA" if lrt is None else repr(lrt), r.get("status", "NA"), dl,
                                    "NA" if diff is None else repr(diff), "NA" if lnl is None else repr(lnl),
                                    r.get("lnl_unconstrained", "NA"), tree])))
open(os.path.join(work, "direct_vs_dndsr.tsv"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
EOF
