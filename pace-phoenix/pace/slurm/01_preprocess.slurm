#!/bin/bash
#SBATCH --job-name=grn-preprocess
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm_%x-%j.out
#SBATCH --error=logs/slurm_%x-%j.err
#
# STEP 1 -- split the raw 10x h5 files into scSAGA-format per-modality files
# (counts.mtx / barcodes.txt / features.txt / pca_50.txt).
#
# Produces:
#   data/3k_rna, data/3k_atac        3k multiome, 2,711 cells each
#   data/10k_rna, data/10k_atac      10k multiome, 11,898 cells each
#   data/6k_rna, data/6k_atac        deterministic 2,711-cell 10k subsample
#                                    (only needed for the 3k+6k workflow)
#   data/atac10k_ext                 external 10k ATAC v1.1, 8,161 cells
#   data/4k_rna                      PBMC 4k v2 (RNA-only baseline)
#
# MEMORY: the 10k multiome h5 is 180,488 features x 11,898 cells (127.5M nonzeros).
# The splitter uses the memory-safe sparse-variance PCA, so 64G is generous.
# TIME: the 10k split took ~7 min on a 10-core laptop; 2 h is ample.
#
#   sbatch pace/slurm/01_preprocess.slurm
#   SKIP_6K=1 SKIP_4K=1 sbatch ...   # skip files your workflow does not need

set -euo pipefail
source "$(dirname "$0")/_common.sh"

cd "$PACE_ROOT"
PY="python"          # the scmint (py3.12) env, activated by _common.sh
S="scripts/00_split_10x_h5.py"

echo "### 3k multiome -> data/3k_{rna,atac}"
"$PY" "$S" --h5 "$PACE_RAW/3k_multiome.h5" --out-subdir 3k

echo "### 10k multiome -> data/10k_{rna,atac}"
"$PY" "$S" --h5 "$PACE_RAW/10k_multiome.h5" --out-subdir 10k

echo "### external 10k ATAC v1.1 -> data/atac10k_ext"
"$PY" "$S" --h5 "$PACE_RAW/atac10k_v1.1.h5" --mode atac --out-subdir atac10k_ext

if [ "${SKIP_6K:-0}" != "1" ]; then
  echo "### 6k: deterministic 2,711-cell subsample of the 10k multiome"
  "$PY" "$S" --h5 "$PACE_RAW/10k_multiome.h5" --out-subdir 6k \
        --subsample 2711 --seed 0
fi

if [ "${SKIP_4K:-0}" != "1" ]; then
  echo "### PBMC 4k v2 -> data/4k_rna"
  "$PY" scripts/00b_prepare_4k_rna.py
fi

echo
echo "=== data layout ==="
for d in "$PACE_DATA"/*/; do
  n=$(ls -1 "$d" | wc -l)
  printf '  %-22s %s files\n' "$(basename "$d")" "$n"
done
du -sh "$PACE_DATA"
