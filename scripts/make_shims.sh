#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-shims}"
TOOLS=(diamond orthofinder mafft muscle blastp blastn makeblastdb)  # add as needed

mkdir -p "${OUTDIR}"
for t in "${TOOLS[@]}"; do
  cat > "${OUTDIR}/${t}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ -n "${RUN_PREFIX:-}" ]]; then
  exec ${RUN_PREFIX} __TOOL__ "$@"
else
  exec __TOOL__ "$@"
fi
EOF
  sed -i.bak "s/__TOOL__/${t}/g" "${OUTDIR}/${t}" && rm -f "${OUTDIR}/${t}.bak"
  chmod +x "${OUTDIR}/${t}"
done
echo "Wrote shims to ${OUTDIR}/"
