#!/bin/bash
#SBATCH --job-name=grn-impute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_%x-%j.out
#SBATCH --error=logs/slurm_%x-%j.err
#
# STEP 3 -- reverse-imputeKNN + the all-cells expression matrix.
#
# For each ATAC query cell: k=20 nearest RNA reference cells in H, weights
# exp(-distance) row-normalised, expression propagated.  Rows of the all-cells
# matrix are real RNA cells then imputed ATAC cells.
#
# MEMORY is the binding constraint, not CPU:
#   unpaired: 20,059 x 36,601 float32  = 2.70 GB  all-cells matrix
#             + 36,601 x 11,898 reference (1.07 GB, float32)
#             + the sparse CPM copy of 10k RNA counts (~0.3 GB)
#   paired  : 23,796 x 36,601 float32  = 3.20 GB
#   A few GB of headroom for scipy's mmread + the dense gene chunks. 96G is
#   deliberately loose; 32G works for the unpaired case but leaves little slack.
#   The dense reference is built in 4,000-gene chunks (GENE_CHUNK) to bound peak
#   memory -- raise it if you have RAM to spare, lower it if the job is OOM-killed.
#
#   sbatch pace/slurm/03_imputation.slurm
#   PACE_EXPERIMENT=paired sbatch pace/slurm/03_imputation.slurm

set -euo pipefail
source "$(dirname "$0")/_common.sh"

cd "$PACE_ROOT"
EXP="${PACE_EXPERIMENT:-unpaired}"

python scripts/05_build_imputation.py --experiment "$EXP"

O="$PACE_ROOT/results/${EXP}_grn"
echo
echo "=== imputation outputs ==="
ls -la "$O"
echo
echo "n_cells: $(cat "$O/n_cells.txt")"
echo "genes  : $(wc -l < "$O/genes.txt")"
