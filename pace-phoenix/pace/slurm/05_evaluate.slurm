#!/bin/bash
#SBATCH --job-name=grn-eval
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm_%x-%j.out
#SBATCH --error=logs/slurm_%x-%j.err
#
# STEP 5 -- evaluate the inferred network against PBMC-TRRUST and PBMC-Blood.
#
# Cheap: reads the network .tsv and scores Precision@K / Recall@K plus the
# deduplicated recovery fractions. 4 CPUs / 16G is generous.
#
#   sbatch pace/slurm/05_evaluate.slurm

set -euo pipefail
source "$(dirname "$0")/_common.sh"

cd "$PACE_ROOT"
EXP="${PACE_EXPERIMENT:-unpaired}"

python scripts/07_evaluate_grn.py --experiment "$EXP"

echo
echo "=== evaluation summary ==="
cat "$PACE_ROOT/results/${EXP}_grn/grn_tfonly_reg/evaluation/evaluation_summary.txt"
