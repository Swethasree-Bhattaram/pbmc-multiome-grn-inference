#!/bin/bash
#SBATCH --job-name=grn-integrate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm_%x-%j.out
#SBATCH --error=logs/slurm_%x-%j.err
#
# STEP 2 -- scSAGA integration: align the ATAC PCA onto the RNA PCA and write the
# joint embedding H that reverse-imputeKNN consumes.
#
#   PACE_EXPERIMENT=unpaired  (default)  rna10k x atac10k_ext   -> H = 20,059 x 30
#   PACE_EXPERIMENT=paired               rna10k x atac10k       -> H = 23,796 x 30
#
# REQUIRES the SAVE-H patch on the scSAGA checkout:
#   python envs/03_patch_scsaga.py        # idempotent; run once after cloning
# The script aborts if joint_embedding_H.npy is not written, so a missing patch
# cannot silently poison the downstream steps.
#
# RESOURCES: the published unpaired run took ~157 s of solver time on a laptop
# (11.8k x 8.2k cells). 8 CPUs / 64G / 2 h is comfortable; scSAGA is CPU-only.
#
#   sbatch pace/slurm/02_integration.slurm
#   PACE_EXPERIMENT=paired sbatch pace/slurm/02_integration.slurm

set -euo pipefail
source "$(dirname "$0")/_common.sh"
export SCAGA_REPO

cd "$PACE_ROOT"
EXP="${PACE_EXPERIMENT:-unpaired}"

echo "### checking the scSAGA SAVE-H patch"
python envs/03_patch_scsaga.py --check

echo "### integration: $EXP"
python scripts/04_run_integration.py --experiment "$EXP"

R="$PACE_ROOT/results/integration_$EXP"
echo
echo "=== integration outputs ==="
ls -la "$R"
echo
echo "=== alignment quality ==="
cat "$R/saga_runtimes.txt"
