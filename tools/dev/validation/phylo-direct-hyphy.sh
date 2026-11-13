#!/usr/bin/env bash
# Runs HyPhy BUSTED directly on the bundled phylogenetic example, without any dndsR code,
# and writes the values dndsR must reproduce. Inputs are prepared here with awk, and the
# trees are written by hand from the design in make-phylo-example.R. JSON is read with
# Python, not dndsR's R parser.
#   tools/dev/validation/phylo-direct-hyphy.sh <hyphy> <example dir> <work dir> [seed] [cpu]
# Writes <work dir>/direct_hyphy.tsv.
set -euo pipefail
hyphy=$1; example=$2; work=$3; seed=${4:-1}; cpu=${5:-1}
mkdir -p "$work"
version=$("$hyphy" --version | grep -o '[0-9]*\.[0-9]*\.[0-9]*' | head -1)

# Hand-written, unrooted, labelled topologies (no branch lengths). Internal nodes are left
# unnamed, so HyPhy names them, unlike dndsR, which names them Node<k>.
declare -A tree=(
  [OG0001]='(((spA{FG},spB{FG}){FG},(spC,spD)),(spE,spF),(spG,spH));'
  [OG0002]='(((spA{FG},spB{FG}){FG},(spC,spD)),(spE,spF),(spG,spH));'
  [OG0003]='(((spA{FG},spB{FG}){FG},(spC,spD)),(spE,spF),spG);'
  [OG0004]='((spA{FG},(spC,spD)),(spE,spF),(spG,spH));'
)

for og in OG0001 OG0002 OG0003 OG0004; do
  d="$work/$og"; mkdir -p "$d"
  # Headers become the taxon; a final TAA, TAG or TGA becomes ---. Sequences are one line each.
  awk '/^>/ { sub(/\|.*/, ""); print; next } { if (match($0, /(TAA|TAG|TGA)$/) && length($0) % 3 == 0) $0 = substr($0, 1, length($0) - 3) "---"; print }' \
    "$example/alignments/$og.fasta" > "$d/alignment.fas"
  printf '%s\n' "${tree[$og]}" > "$d/tree.nwk"
  ( cd "$d" && "$hyphy" CPU="$cpu" ENV="RANDOM_SEED=$seed" busted --alignment alignment.fas --tree tree.nwk \
      --branches FG --code Universal --srv Yes --rates 3 --syn-rates 3 --multiple-hits None --error-sink No \
      --grid-size 250 --starting-points 1 --output busted.json > hyphy.log 2>&1 ) &
done
wait

python3 - "$work" "$version" "$seed" "$cpu" <<'EOF'
import json, sys, os
work, version, seed, cpu = sys.argv[1:]
rows = ["orthogroup\thyphy_version\tseed\tcpu\tlrt\tp_value\tlnl_unconstrained\tlnl_constrained\tn_tested_branches"]
for og in ["OG0001", "OG0002", "OG0003", "OG0004"]:
    j = json.load(open(os.path.join(work, og, "busted.json")))
    fits = j["fits"]
    con = fits.get("Constrained model", {}).get("Log Likelihood")
    tested = sum(1 for v in j["tested"]["0"].values() if v == "test")
    values = [j["test results"]["LRT"], j["test results"]["p-value"], fits["Unconstrained model"]["Log Likelihood"], con]
    rows.append("\t".join([og, version, seed, cpu] + ["NA" if v is None else repr(v) for v in values] + [str(tested)]))
open(os.path.join(work, "direct_hyphy.tsv"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
EOF
