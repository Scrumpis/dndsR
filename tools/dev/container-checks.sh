#!/usr/bin/env bash
# Validates a built dndsR image from the inside: tools, documented CLI and R workflows,
# failure behavior, and expected example results. The same script runs under Docker and
# Apptainer. Run it in an empty writable folder, as the current user:
#   docker run --rm -u "$(id -u):$(id -g)" -e HOME=/work -v "$PWD/work:/work" -w /work \
#     -v "$PWD/tools/dev/container-checks.sh:/checks.sh:ro" IMAGE bash /checks.sh [expected version]
set -u
expected_version=${1:-}
failures=0
log=$(mktemp)

# In GitHub Actions, failures also appear as public annotations, since logs need a login.
report() {
  echo "FAIL  $1"
  if [ -n "${GITHUB_ACTIONS:-}" ]; then
    printf '::error title=container check failed::%s%%0A%s\n' "$1" \
      "$(tail -n 12 "$log" | sed 's/%/%25/g' | tr '\n' '~' | sed 's/~/%0A/g' | cut -c1-1500)"
  fi
  failures=$((failures + 1))
}
pass() { # description, command...
  local d=$1; shift
  if "$@" >"$log" 2>&1; then echo "ok    $d"; else report "$d"; tail -n 12 "$log"; fi
}
fail_expected() {
  local d=$1; shift
  if "$@" >"$log" 2>&1; then report "$d (should have failed)"; else echo "ok    $d"; fi
}
has_text() { grep -q -- "$1" "$2"; }

echo "== Tools and dependencies"
pass "dndsr is on PATH" command -v dndsr
pass "DIAMOND runs" diamond version
pass "OrthoFinder runs" orthofinder -h
pass "required R packages load" Rscript -e 'for (p in c("dndsR","orthologr","Biostrings","GenomicFeatures","txdbmaker","pwalign","topGO","GO.db","IHW","qvalue","RIdeogram","lme4","testthat")) stopifnot(requireNamespace(p, quietly = TRUE))'
if [ -n "$expected_version" ]; then
  pass "dndsR version is $expected_version" Rscript -e "stopifnot(as.character(packageVersion('dndsR')) == '$expected_version')"
fi
Rscript -e 'ip <- installed.packages()[, c("Package", "Version")]; print(ip[ip[, "Package"] %in% c("dndsR","orthologr","Biostrings","GenomicFeatures"), ], quote = FALSE); cat(R.version.string, "\n")'

echo "== Documented CLI workflow"
pass "help lists commands" sh -c 'dndsr --help | grep -q verify_example'
pass "command help works" sh -c 'dndsr calculate_dnds --help | grep -q comparison_file'
pass "copy the example" dndsr example -o cli
cd cli || exit 1
pass "extract_cds" dndsr extract_cds -c comparison_file.txt
pass "calculate_dnds" dndsr calculate_dnds -c comparison_file.txt -t 2
pass "verify_example: 14 checks" sh -c 'dndsr verify_example 2>&1 | grep -q "All 14 example checks passed"'
pass "append_annotations" dndsr append_annotations -c comparison_file.txt
pass "term_enrichment" dndsr term_enrichment -c comparison_file.txt
pass "term output" test -s example_AvB/example_AvB_q_KEGG_enrichment.tsv
pass "ipr_enrichment" dndsr ipr_enrichment -c comparison_file.txt --max-prop 1
pass "InterPro output" test -s example_AvB/example_AvB_q_IPR_Domain_enrichment.tsv
pass "go_enrichment" dndsr go_enrichment -c comparison_file.txt --node-min 1
pass "GO output" sh -c 'ls example_AvB/*_q_GO_BP_*.tsv >/dev/null'
pass "dnds_ideogram" dndsr dnds_ideogram -c comparison_file.txt
pass "ideogram output" sh -c 'ls example_AvB/*ideogram.svg >/dev/null'
pass "dnds_summary with regions" dndsr dnds_summary -c comparison_file.txt --regions-bed regions.bed
pass "two comparisons" sh -c 'dndsr extract_cds -c two_comparisons.txt && dndsr calculate_dnds -c two_comparisons.txt -t 2 && dndsr append_annotations -c two_comparisons.txt'
pass "dnds_contrast" dndsr dnds_contrast -c two_comparisons.txt --contrast-file contrast_file.txt
pass "contrast output" test -s contrast/AvB_vs_AvC_global_q_vs_q_dnds_contrast.tsv
cd .. || exit 1

echo "== Failure behavior"
fail_expected "unknown command fails" dndsr nonsense
fail_expected "unknown option fails" dndsr extract_cds -c x --overwite
fail_expected "missing input fails" dndsr calculate_dnds -c does-not-exist.txt
fail_expected "gene-level state contrast fails" dndsr dnds_state_contrast --level gene --dnds-annot-files cli/example_AvB/example_AvB_dnds_annot.tsv

echo "== Documented R workflow"
cat >r-workflow.R <<'RS'
library(dndsR)
folder <- dndsr_example("r_example")
setwd(folder)
extract_cds(comparison_file = "comparison_file.txt")
calculate_dnds(comparison_file = "comparison_file.txt", comp_cores = 2)
checks <- verify_example()
stopifnot(all(checks$passed), nrow(checks) == 14L)
append_annotations(comparison_file = "comparison_file.txt")
dnds <- read.delim("example_AvB/example_AvB_dnds.tsv")
stopifnot(nrow(dnds) == 32L)
RS
pass "R functions, end to end" Rscript r-workflow.R

echo
if [ "$failures" -gt 0 ]; then
  echo "$failures container check(s) failed"
  exit 1
fi
echo "All container checks passed"
