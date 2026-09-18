#!/bin/bash
# Build the per-dataset files the pipeline needs, INCLUDING the PCA.
#
#   bash scripts/prepare_data.sh            # everything
#   bash scripts/prepare_data.sh 3k 10k     # only some datasets
#
# This is the PCA step: for each dataset it writes
#     data/<name>/pca_50.txt    cells x 50 PCs, on log1p(CPM/1e4) of the
#                               top-2000 variable features
# plus counts.mtx / barcodes.txt / features.txt where an RNA reference is needed.
#
# scSAGA consumes ONLY pca_50.txt per dataset (it reads the `pca` key and nothing
# else), so if you already have PCs for a dataset you can skip this entirely and
# point config.yml straight at those files:
#
#     my_rna:
#       modality: rna
#       dir:  /path/to/expr          # counts.mtx barcodes.txt features.txt
#       pca:  /path/to/my_pcs.txt    # any cells x N file, np.loadtxt-readable
#
# Run in the scmint env (needs h5py, sklearn).
set -euo pipefail

ROOT="${PIPELINE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RAW="${PACE_RAW:-$ROOT/raw}"
cd "$ROOT"
mkdir -p logs

# which datasets to build (default: all that this repo knows how to build)
WANT=("$@")
want() {
  [ ${#WANT[@]} -eq 0 ] && return 0
  for w in "${WANT[@]}"; do [ "$w" = "$1" ] && return 0; done
  return 1
}

log() { echo; echo "### $*"; }

# --------------------------------------------------------------------------- #
# 3k and 10k multiome: one h5 each, both modalities out
# --------------------------------------------------------------------------- #
for tag in 3k 10k; do
  want "$tag" || continue
  H5="$RAW/${tag}_multiome.h5"
  if [ ! -f "$H5" ]; then
    echo "SKIP $tag: no $H5"; continue
  fi
  log "split $tag multiome (+ PCA for both modalities)"
  python scripts/split_10x_h5.py --h5 "$H5" --out-subdir "$tag" --data-dir "$ROOT/data"
done

# --------------------------------------------------------------------------- #
# external 10k ATAC v1.1 (ATAC only)
# --------------------------------------------------------------------------- #
if want atac10k_ext && [ -f "$RAW/atac10k_v1.1.h5" ]; then
  log "split external 10k ATAC v1.1 (+ PCA)"
  python scripts/split_10x_h5.py --h5 "$RAW/atac10k_v1.1.h5" --mode atac \
        --out-subdir atac10k_ext --data-dir "$ROOT/data"
fi

# --------------------------------------------------------------------------- #
# 6k = deterministic 2,711-cell subsample of the 10k multiome (same cells both
# modalities).  Only needed by the `ab_6k` experiment in config.yml.
# --------------------------------------------------------------------------- #
if want 6k && [ -f "$RAW/10k_multiome.h5" ]; then
  log "split 6k (2,711-cell subsample of 10k)"
  python scripts/split_10x_h5.py --h5 "$RAW/10k_multiome.h5" --out-subdir 6k \
        --subsample 2711 --seed 0 --data-dir "$ROOT/data"
fi

# --------------------------------------------------------------------------- #
# PBMC 4k v2 tar.gz -> data/4k_rna (RNA only; used by the RNA-only baseline)
# --------------------------------------------------------------------------- #
if want 4k && [ -f "$RAW/pbmc4k_filtered_gene_bc_matrices.tar.gz" ]; then
  log "prepare PBMC 4k v2 -> data/4k_rna"
  python scripts/prepare_4k_rna.py
fi

# --------------------------------------------------------------------------- #
log "dataset inventory"
for d in "$ROOT"/data/*/; do
  n=$(basename "$d")
  pca="-"; [ -f "$d/pca_50.txt" ] && pca="pca_50.txt"
  cnt="-"; [ -f "$d/counts.mtx" ] && cnt="counts.mtx"
  cells="-"
  [ -f "$d/barcodes.txt" ] && cells="$(wc -l < "$d/barcodes.txt") cells"
  printf '  %-14s %-12s %-12s %s\n' "$n" "$pca" "$cnt" "$cells"
done
