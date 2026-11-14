#!/usr/bin/env bash
# Freezes stage 2 (docs/dev/phylo-calibration-stage2.md): the committed code, a dedicated
# read-only dndsR library installed from it, and the image by digest.
#   tools/dev/validation/phylo-stage2-freeze.sh <results root> <image> <expected sha256> <params.json> <fg_segments.tsv>
# Run from a clean checkout. Refuses if the tree has uncommitted changes, or if the root is
# already frozen at another commit.
set -euo pipefail
root=$(realpath -m "$1"); sif=$(realpath "$2"); want=$3; params=$(realpath "$4"); segs=$(realpath "$5")
repo=$(git rev-parse --show-toplevel)
[ -z "$(git -C "$repo" status --porcelain --untracked-files=no)" ] || { echo "uncommitted changes; commit first" >&2; exit 1; }
commit=$(git -C "$repo" rev-parse HEAD)
if [ -f "$root/freeze.json" ]; then
  grep -q "\"commit\": \"$commit\"" "$root/freeze.json" && { echo "already frozen at $commit"; exit 0; }
  echo "$root is frozen at another commit; use a new root" >&2; exit 1
fi
got=$(sha256sum "$sif" | cut -d' ' -f1)
[ "$got" = "$want" ] || { echo "image sha256 $got, expected $want" >&2; exit 1; }
mkdir -p "$root/src" "$root/lib" "$root/inputs"
git -C "$repo" archive "$commit" | tar -x -C "$root/src"
cp "$params" "$root/inputs/params.json"; cp "$segs" "$root/inputs/fg_segments.tsv"
singularity exec --cleanenv --bind "$root" "$sif" R CMD INSTALL --no-test-load -l "$root/lib" "$root/src" > "$root/lib_install.log" 2>&1
chmod -R a-w "$root/lib" "$root/src"
python3 -I - "$root" "$commit" "$sif" "$got" <<'EOF'
import hashlib, json, os, sys, datetime
root, commit, sif, sha = sys.argv[1:]
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
json.dump({"commit": commit, "image": sif, "image_sha256": sha, "library": os.path.join(root, "lib"),
           "source": os.path.join(root, "src"), "params_md5": md5(os.path.join(root, "inputs/params.json")),
           "fg_segments_md5": md5(os.path.join(root, "inputs/fg_segments.tsv")),
           "frozen_utc": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")},
          open(os.path.join(root, "freeze.json"), "w"), indent=1)
EOF
cat "$root/freeze.json"
