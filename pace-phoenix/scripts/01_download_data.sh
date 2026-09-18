#!/bin/bash
# Download the raw 10x datasets used by the scSAGA + reverse-imputeKNN + Arboreto
# GRN pipeline, and verify each file's MD5 before it is used.
#
# Verified against the copies used for the published runs in this repo
# (checksums recorded 2026-09-18; all four files re-downloaded and hashed).
#
# Usage:
#   bash 01_download_data.sh                 # into $PACE_RAW (default <workspace>/raw)
#   PACE_RAW=/scratch/$USER/grn/raw bash 01_download_data.sh
#   FORCE=1 bash 01_download_data.sh         # re-download even if present + valid
#
# Run this on a LOGIN node (it is I/O only, ~316 MB), or as a batch job with
# `sbatch pace/slurm/00_download.slurm`.
set -euo pipefail

PACE_ROOT="${PACE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RAW="${PACE_RAW:-$PACE_ROOT/raw}"
FORCE="${FORCE:-0}"
mkdir -p "$RAW"
cd "$RAW"

BASE="https://cf.10xgenomics.com/samples"

# name|url|md5|bytes
FILES=(
  "3k_multiome.h5|${BASE}/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5|e326066b51ec8975197c29a7f911a4fd|38844318"
  "10k_multiome.h5|${BASE}/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5|df86844b99161b9487090d91e644745e|192125528"
  "atac10k_v1.1.h5|${BASE}/cell-atac/1.1.0/atac_pbmc_10k_v1/atac_pbmc_10k_v1_filtered_peak_bc_matrix.h5|5ee75b4d7b5d70945a3c33a78c63582d|67041999"
  "pbmc4k_filtered_gene_bc_matrices.tar.gz|${BASE}/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz|f61f4deca423ef0fa82d63fdfa0497f7|18423814"
)

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}';
  else md5 -q "$1"; fi
}

fail=0
for row in "${FILES[@]}"; do
  IFS='|' read -r name url want_md5 want_bytes <<< "$row"
  echo "=== $name"

  if [ -f "$name" ] && [ "$FORCE" != "1" ]; then
    have=$(md5_of "$name")
    if [ "$have" = "$want_md5" ]; then
      echo "    already present and verified (md5 $have)"
      continue
    fi
    echo "    present but md5 mismatch (have $have, want $want_md5) -> re-downloading"
  fi

  # -C - resumes a partial download from an earlier attempt
  curl -fSL --retry 5 --retry-delay 5 -C - -o "$name" "$url"
  have_bytes=$(wc -c < "$name" | tr -d ' ')
  have=$(md5_of "$name")

  if [ "$have" != "$want_md5" ]; then
    echo "    FAIL md5: have $have want $want_md5" >&2
    fail=1
    continue
  fi
  if [ "$have_bytes" != "$want_bytes" ]; then
    echo "    WARN size: have $have_bytes want $want_bytes" >&2
  fi
  echo "    OK  $have_bytes bytes  md5 $have"
done

if [ "$fail" != "0" ]; then
  echo "ERROR: one or more downloads failed verification" >&2
  exit 1
fi

cat > "$RAW/CHECKSUMS.txt" <<'EOF'
MD5 of the raw 10x inputs used by the pbmc-multiome-grn-inference pipeline.
Verified 2026-09-18.  See 01_download_data.sh for the source URLs.

e326066b51ec8975197c29a7f911a4fd  3k_multiome.h5
df86844b99161b9487090d91e644745e  10k_multiome.h5
5ee75b4d7b5d70945a3c33a78c63582d  atac10k_v1.1.h5
f61f4deca423ef0fa82d63fdfa0497f7  pbmc4k_filtered_gene_bc_matrices.tar.gz
EOF

echo
echo "All raw inputs verified in $RAW"
ls -la
