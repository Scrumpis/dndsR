#!/usr/bin/env bash
# Inputs for running codeml, aBSREL, and RELAX directly on the bundled alignment example,
# without dndsR code: PHYLIP files from awk, and trees written by hand from the design in
# make-phylo-example.R. Writes into <example>/truth/direct/<OG>/, beside the BUSTED inputs.
# The tests (test-phylo-lineage.R) run the engines on these files and compare with dndsR.
#   tools/dev/validation/phylo-direct-lineage.sh <example dir>
set -euo pipefail
example=$1
direct="$example/truth/direct"

# Unrooted topologies. OG0001 has every taxon; OG0004 lacks spB, so the foreground clade
# (spA, spB and their stem) collapses onto spA's branch.
declare -A plain=(
  [OG0001]='(((spA,spB),(spC,spD)),(spE,spF),(spG,spH));'
  [OG0004]='((spA,(spC,spD)),(spE,spF),(spG,spH));'
)
declare -A codeml_fg=(
  [OG0001]='(((spA #1,spB #1) #1,(spC,spD)),(spE,spF),(spG,spH));'
  [OG0004]='((spA #1,(spC,spD)),(spE,spF),(spG,spH));'
)
declare -A absrel_all=(
  [OG0001]='(((spA{FG},spB{FG}){FG},(spC{FG},spD{FG}){FG}){FG},(spE{FG},spF{FG}){FG},(spG{FG},spH{FG}){FG});'
  [OG0004]='((spA{FG},(spC{FG},spD{FG}){FG}){FG},(spE{FG},spF{FG}){FG},(spG{FG},spH{FG}){FG});'
)
declare -A relax=(
  [OG0001]='(((spA{TEST},spB{TEST}){TEST},(spC{REFERENCE},spD{REFERENCE}){REFERENCE}){REFERENCE},(spE{REFERENCE},spF{REFERENCE}){REFERENCE},(spG{REFERENCE},spH{REFERENCE}){REFERENCE});'
  [OG0004]='((spA{TEST},(spC{REFERENCE},spD{REFERENCE}){REFERENCE}){REFERENCE},(spE{REFERENCE},spF{REFERENCE}){REFERENCE},(spG{REFERENCE},spH{REFERENCE}){REFERENCE});'
)

for og in OG0001 OG0004; do
  d="$direct/$og"
  n=$(grep -c '^>' "$d/alignment.fas")
  # PHYLIP for codeml from the stored FASTA (already one line per sequence, stops masked).
  awk -v n="$n" 'BEGIN { getline; } /^>/ { name = substr($0, 2); next } { seqs[++k] = name "  " $0; len = length($0) }
       END { printf " %d %d\n", n, len; for (i = 1; i <= k; i++) print seqs[i] }' <(echo; cat "$d/alignment.fas") > "$d/alignment.phy"
  printf ' %d 1\n%s\n' "$n" "${plain[$og]}" > "$d/tree_codeml.nwk"
  printf ' %d 1\n%s\n' "$n" "${codeml_fg[$og]}" > "$d/tree_codeml_fg.nwk"
  printf '%s\n' "${absrel_all[$og]}" > "$d/tree_absrel.nwk"
  printf '%s\n' "${relax[$og]}" > "$d/tree_relax.nwk"
done
echo "Wrote direct inputs for OG0001 and OG0004 in $direct"
