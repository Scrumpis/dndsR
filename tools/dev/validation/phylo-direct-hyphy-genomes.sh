#!/usr/bin/env bash
# Runs HyPhy BUSTED directly on orthogroups of the synthetic genome example, without dndsR:
# alignments come from the designed coding sequences (truth/designed_cds.fasta), which have
# no indels in these orthogroups, with the final stop removed by awk; trees are written by
# hand. JSON is read with Python.
#   tools/dev/validation/phylo-direct-hyphy-genomes.sh <hyphy> <example dir> <work dir> [seed] [cpu]
# Writes <work dir>/direct_hyphy.tsv: OG01 and OG05 on the species tree, and OG04 and OG05 on
# their gene trees.
set -euo pipefail
hyphy=$1; example=$2; work=$3; seed=${4:-1}; cpu=${5:-1}
mkdir -p "$work"
version=$("$hyphy" --version | grep -o '[0-9]*\.[0-9]*\.[0-9]*' | head -1)

# run name, orthogroup, hand-written labelled topology
runs=(
  "OG01_species OG01 ((((spA{FG},spB{FG}){FG},spC),(spD,spE)),spF);"
  "OG05_species OG05 ((((spA{FG},spB{FG}){FG},spC),(spD,spE)),spF);"
  "OG04_gene OG04 ((((spA{FG},spB{FG}){FG},spC),(spD,spF)),spE);"
  "OG05_gene OG05 (((spA{FG},spB{FG}){FG},(spC,spD)),spE,spF);"
)
for r in "${runs[@]}"; do
  read -r name og tree <<<"$r"
  d="$work/$name"; mkdir -p "$d"
  # Genes of the orthogroup, from the long table; sequences from the designed CDS, final stop removed.
  awk -F'\t' -v og="$og" 'NR > 1 && $1 == og { print $3 "\t" $2 }' "$example/orthogroup_table.tsv" > "$d/genes.tsv"
  awk -F'\t' 'FNR == NR { taxon[$1] = $2; next }
              /^>/ { id = substr($1, 2); keep = (id in taxon); if (keep) print ">" taxon[id]; next }
              keep { print substr($0, 1, length($0) - 3) }' "$d/genes.tsv" "$example/truth/designed_cds.fasta" > "$d/alignment.fas"
  printf '%s\n' "$tree" > "$d/tree.nwk"
  ( cd "$d" && "$hyphy" CPU="$cpu" ENV="RANDOM_SEED=$seed" busted --alignment alignment.fas --tree tree.nwk \
      --branches FG --code Universal --srv Yes --rates 3 --syn-rates 3 --multiple-hits None --error-sink No \
      --grid-size 250 --starting-points 1 --output busted.json > hyphy.log 2>&1 ) &
done
wait

python3 - "$work" "$version" "$seed" "$cpu" <<'EOF'
import json, os, sys
work, version, seed, cpu = sys.argv[1:]
rows = ["run\torthogroup\ttopology\thyphy_version\tseed\tcpu\tlrt\tp_value\tlnl_unconstrained\tlnl_constrained\ttested_branches"]
for name in ["OG01_species", "OG05_species", "OG04_gene", "OG05_gene"]:
    j = json.load(open(os.path.join(work, name, "busted.json")))
    fits = j["fits"]
    con = fits.get("Constrained model", {}).get("Log Likelihood")
    tested = sorted(k for k, v in j["tested"]["0"].items() if v == "test")
    values = [j["test results"]["LRT"], j["test results"]["p-value"], fits["Unconstrained model"]["Log Likelihood"], con]
    og, topo = name.split("_")
    rows.append("\t".join([name, og, topo, version, seed, cpu] + ["NA" if v is None else repr(v) for v in values] + [str(len(tested))]))
open(os.path.join(work, "direct_hyphy.tsv"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
EOF
