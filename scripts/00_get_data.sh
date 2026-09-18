#!/bin/bash
# =============================================================================
#  DOWNLOAD + FORMAT the 10x datasets  ->  data/<name>/{pca_50.txt, counts.mtx,
#                                                barcodes.txt, features.txt}
# =============================================================================
#  One script, start to finish.  Run it once on a login node.
#
#      bash scripts/00_get_data.sh                 # download + format everything
#      bash scripts/00_get_data.sh download         # only fetch the raw files
#      bash scripts/00_get_data.sh format           # only build data/<name>/
#      FORCE=1 bash scripts/00_get_data.sh          # re-download + rebuild
#
#  WHAT YOU GET -- 4 raw files become 6 formatted datasets:
#
#      raw file                            ->  data/<name>/            modality
#      ------------------------------------------------------------------
#      pbmc_granulocyte_sorted_3k  .h5     ->  3k_rna, 3k_atac        both
#      pbmc_granulocyte_sorted_10k .h5     ->  10k_rna, 10k_atac      both
#      atac_pbmc_10k_v1            .h5     ->  atac10k_ext           ATAC
#      pbmc4k_..._matrices.tar.gz          ->  4k_rna                RNA
#      (10k multiome, 2,711-cell subsample)->  6k_rna, 6k_atac        both
#
#  That is the 10k multiome, the 10k ATAC-seq PBMC, plus the 3k multiome and the
#  4k RNA-seq dataset, and the 6k subsample used by the `ab_6k` experiment.
#
#  Per dataset, `pca_50.txt` is 50 PCs on log1p(CPM/1e4) of the top-2000 variable
#  features.  scSAGA consumes ONLY that file.  counts.mtx / barcodes.txt /
#  features.txt are written too, but only matter when the dataset is used as an
#  RNA REFERENCE (i.e. its real expression is propagated onto ATAC cells).
#
#  ALREADY HAVE THE DATA?  Skip this script completely and point config.yml at
#  what you have -- any cells x N PC file works:
#
#      my_rna:
#        modality: rna
#        dir:  /path/to/expr          # counts.mtx barcodes.txt features.txt
#        pca:  /path/to/my_pcs.txt    # your own PCs
#
#  Needs the scmint env (h5py, sklearn, scipy).
# =============================================================================
set -euo pipefail

ROOT="${PIPELINE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RAW="${PACE_RAW:-$ROOT/raw}"
FORCE="${FORCE:-0}"
STAGE="${1:-all}"
cd "$ROOT"
mkdir -p "$RAW" logs

BASE="https://cf.10xgenomics.com/samples"
# name|url|md5|bytes   (checksums verified 2026-09-18 against the published runs)
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

# --------------------------------------------------------------------------- #
# 1. download (idempotent: skips files already present with a matching md5)
# --------------------------------------------------------------------------- #
do_download() {
  echo "### download raw 10x inputs -> $RAW"
  local fail=0
  for row in "${FILES[@]}"; do
    IFS='|' read -r name url want_md5 want_bytes <<< "$row"
    echo "=== $name"

    if [ -f "$RAW/$name" ] && [ "$FORCE" != "1" ]; then
      have=$(md5_of "$RAW/$name")
      if [ "$have" = "$want_md5" ]; then
        echo "    present + verified ($have)"; continue
      fi
      echo "    md5 mismatch (have $have want $want_md5) -> re-downloading"
    fi

    # -C - resumes an interrupted download instead of starting over
    curl -fSL --retry 5 --retry-delay 5 -C - -o "$RAW/$name" "$url"

    have=$(md5_of "$RAW/$name")
    if [ "$have" != "$want_md5" ]; then
      echo "    FAIL md5: have $have want $want_md5" >&2; fail=1; continue
    fi
    echo "    OK  $(wc -c < "$RAW/$name" | tr -d ' ') bytes  md5 $have"
  done

  [ "$fail" = "0" ] || { echo "ERROR: verification failed" >&2; exit 1; }
  cat > "$RAW/CHECKSUMS.txt" <<'EOF'
MD5 of the raw 10x inputs used by the pbmc-multiome-grn-inference pipeline.
Verified 2026-09-18.  Source URLs are in scripts/00_get_data.sh.

e326066b51ec8975197c29a7f911a4fd  3k_multiome.h5
df86844b99161b9487090d91e644745e  10k_multiome.h5
5ee75b4d7b5d70945a3c33a78c63582d  atac10k_v1.1.h5
f61f4deca423ef0fa82d63fdfa0497f7  pbmc4k_filtered_gene_bc_matrices.tar.gz
EOF
  echo "    all raw inputs verified in $RAW"
}

# --------------------------------------------------------------------------- #
# 2. format  (raw files -> data/<name>/ with pca_50.txt + counts)
# --------------------------------------------------------------------------- #
do_format() {
  echo "### format -> data/<name>/"

  # one h5, both modalities
  for tag in 3k 10k; do
    [ -f "$RAW/${tag}_multiome.h5" ] || { echo "SKIP $tag: no raw h5"; continue; }
    echo "=== $tag multiome -> ${tag}_rna + ${tag}_atac (+ PCA)"
    python scripts/split_10x_h5.py --h5 "$RAW/${tag}_multiome.h5" \
           --out-subdir "$tag" --data-dir "$ROOT/data"
  done

  # external 10k ATAC-seq (ATAC only)
  if [ -f "$RAW/atac10k_v1.1.h5" ]; then
    echo "=== external 10k ATAC v1.1 -> atac10k_ext (+ PCA)"
    python scripts/split_10x_h5.py --h5 "$RAW/atac10k_v1.1.h5" --mode atac \
           --out-subdir atac10k_ext --data-dir "$ROOT/data"
  fi

  # 6k = deterministic 2,711-cell subsample of the 10k multiome, same cells in
  # both modalities.  Only used by the `ab_6k` experiment.
  if [ -f "$RAW/10k_multiome.h5" ]; then
    echo "=== 6k subsample (2,711 cells of 10k, seed 0) -> 6k_rna + 6k_atac"
    python scripts/split_10x_h5.py --h5 "$RAW/10k_multiome.h5" --out-subdir 6k \
           --subsample 2711 --seed 0 --data-dir "$ROOT/data"
  fi

  # PBMC 4k v2 tar.gz -> RNA only (RNA-only baseline; no PCA)
  if [ -f "$RAW/pbmc4k_filtered_gene_bc_matrices.tar.gz" ]; then
    echo "=== PBMC 4k v2 -> data/4k_rna"
    python scripts/prepare_4k_rna.py
  fi

  echo
  echo "### inventory"
  printf '  %-14s %-12s %-12s %s\n' DATASET PCA COUNTS CELLS
  for d in "$ROOT"/data/*/; do
    [ -d "$d" ] || continue
    n=$(basename "$d"); pca="-"; cnt="-"; cells="-"
    [ -f "$d/pca_50.txt" ] && pca="pca_50.txt"
    [ -f "$d/counts.mtx" ] && cnt="counts.mtx"
    [ -f "$d/barcodes.txt" ] && cells="$(wc -l < "$d/barcodes.txt") cells"
    printf '  %-14s %-12s %-12s %s\n' "$n" "$pca" "$cnt" "$cells"
  done
  echo
  echo "note: a dataset with no pca_50.txt cannot enter the integration."
  echo "      A dataset with no counts.mtx cannot serve as an RNA reference."
}

case "$STAGE" in
  all)      do_download; echo; do_format ;;
  download) do_download ;;
  format)   do_format ;;
  *) echo "usage: $0 [all|download|format]" >&2; exit 2 ;;
esac

echo
echo "done. next:  python pipeline/run.py --list   then   --check"
