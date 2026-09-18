#!/bin/bash
#SBATCH --job-name=grn-grnboost2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/slurm_%x-%j.out
#SBATCH --error=logs/slurm_%x-%j.err
#
# STEP 4 -- Arboreto GRNBoost2 (the expensive step).
#
# Input is the all-cells matrix subset to the UNION of trrust_tf.txt target genes
# and the tf_only.txt regulators: 20,059 cells x 2,852 columns for the unpaired
# run (~0.23 GB float32), which is why 64G is plenty here.
#
# PARALLELISM -- yes, this scales with workers.  Arboreto accepts a pre-built
# dask Client (client_or_address), so the worker count is a runtime knob.
# Measured on the real matrix shape (2,000 cells x 400 genes, seed 666; identical
# 21,998 edges at every setting):
#     1 worker   92.0 s
#     2 workers  47.2 s
#     4 workers  25.2 s
#     8 workers  18.6 s
# Keep GRN_WORKERS * GRN_THREADS_PER_WORKER == --cpus-per-task.  With >8 workers
# the gain flattens: work is split per TARGET GENE, so the slowest gene bounds it.
# The published unpaired run (20,059 x 2,852, 4 workers) took ~2 h 10 m.
#
# MULTI-NODE: set GRN_JOBQUEUE=1 to let this job submit its own worker jobs to
# SLURM (dask_jobqueue.SLURMCluster), e.g. 2 nodes x 8 cores:
#     GRN_JOBQUEUE=1 GRN_NODE_WORKERS=2 sbatch pace/slurm/04_grn.slurm
# That requires dask-jobqueue in the grn39 env (installed by 02_setup_envs.sh)
# and the submitting node being allowed to call sbatch.
#
#   sbatch pace/slurm/04_grn.slurm
#   GRN_WORKERS=4 sbatch pace/slurm/04_grn.slurm
#   GRN_THREADS_PER_WORKER=2 GRN_WORKERS=4 sbatch pace/slurm/04_grn.slurm

set -euo pipefail
ENV_NAME="${GRN_ENV:-grn39}"          # Arboreto needs the python 3.9 env
source "$(dirname "$0")/_common.sh"

cd "$PACE_ROOT"
EXP="${PACE_EXPERIMENT:-unpaired}"

# default: one worker process per core
export GRN_WORKERS="${GRN_WORKERS:-$SLURM_CPUS_PER_TASK}"
export GRN_CPUS="${GRN_CPUS:-$SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$GRN_THREADS_PER_WORKER}"

echo "### GRN: experiment=$EXP workers=$GRN_WORKERS " \
     "threads_per_worker=${GRN_THREADS_PER_WORKER:-1} jobqueue=${GRN_JOBQUEUE:-0}"

python scripts/06_run_grn.py --experiment "$EXP"

O="$PACE_ROOT/results/${EXP}_grn/grn_tfonly_reg"
echo
echo "=== run info ==="
cat "$O/run_info.txt"
echo
echo "=== outputs ==="
ls -la "$O"
