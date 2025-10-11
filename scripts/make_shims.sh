#!/usr/bin/env bash
set -euo pipefail

SIF="${1:-dndsr.sif}"
OUTDIR="${2:-shims}"

TOOLS=(
  diamond
  orthofinder
  mafft
  muscle
  blastp
  blastn
  makeblastdb
  # add more if your pipeline calls them
)

mkdir -p "${OUTDIR}"

ENGINE=""
if command -v singularity >/dev/null 2>&1; then ENGINE="singularity"
elif command -v apptainer  >/dev/null 2>&1; then ENGINE="apptainer"
else
  echo "No singularity/apptainer in PATH" >&2
  exit 1
fi

for t in "${TOOLS[@]}"; do
  cat > "${OUTDIR}/${t}" <<EOF
#!/usr/bin/env bash
exec ${ENGINE} exec -B "\$PWD":/work "${SIF}" ${t} "\$@"
EOF
  chmod +x "${OUTDIR}/${t}"
done

echo "Wrote shims to ${OUTDIR}/"
