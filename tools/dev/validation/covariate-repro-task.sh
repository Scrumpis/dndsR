#!/bin/bash
# One run of the reproduction: the same labels through the release candidate and through covariate-refit.py,
# at pos_threshold 1 and q0.95, beside the terms from the original confirmation run.
#   covariate-repro-task.sh <manifest> <line>
# Needs RT (runtime folder), V (folder with lib and go_ancestors.tsv), SIF, CODE (tools/dev/validation), PY (python with numpy).
read -r pair sc seed <<< "$(sed -n "${2}p" "$1")"
case $pair in
  AtAl) A=$RT/final/new/AtAl ;;
  SlSt) A=$RT/v1/pairs/SlSt/SlSt ;;
esac
stored=$RT/v1/confirm/$pair/$sc
O=$V/repro/$pair/$sc; mkdir -p "$O"
case $sc in
  permute|rotate|strata5|strata10) L=$sc:$seed ;;
  *) L=$stored/sim_$seed.tsv.gz ;;
esac
SI="singularity exec --cleanenv --bind $RT --bind $V --bind $CODE --env R_LIBS_USER=$V/lib --env OMP_NUM_THREADS=1"
for t in 1 q0.95; do
  tab=$O/tab_${t}_$seed.tsv
  $SI --env DNDSR_POS_THRESHOLD=$t --env DNDSR_KEEP_TABLE="$tab" "$SIF" Rscript --vanilla "$CODE/covariate-confirm.R" "$A/${pair}_dnds_annot.tsv" "$L" "$O/terms_${t}_$seed.tsv" >/dev/null 2>&1
  "$PY" -I "$CODE/covariate-refit.py" "$tab" "$V/go_ancestors.tsv" "$O/py_${t}_$seed.tsv" "$t" || echo "refit failed $pair $sc $seed $t"
  rm -f "$tab"
done
echo "done $pair $sc $seed"
