#!/bin/bash
#SBATCH --job-name=grn-grnboost2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_%x-%j.out
#SBATCH --error=logs/slurm_%x-%j.err
#
# STEP 4 -- Arboreto GRNBoost2 (the expensive step).
#
# ---------------------------------------------------------------------------
# WORKERS: set to 16 by default (see GRN_WORKERS below).  --cpus-per-task=16 and
# GRN_WORKERS=16 must MATCH: LocalCluster(n_workers=16, threads_per_worker=1)
# needs exactly 16 cores.  If you change one, change the other.
#
# MEMORY: 64G is sized for 16 workers and is generous.  Per-worker cost is
# dominated by the TF matrix that Arboreto broadcasts to EVERY worker
# (client.scatter(broadcast=True) in arboreto/core.py create_graph):
#
#     n_TFs x n_cells x 4 bytes
#     unpaired: 816 TFs x 20,059 cells = 65 MB per worker
#     with 16 workers: 16 x 65 MB = ~1.0 GB total for the broadcast copy
#
#   plus per-worker interpreter + numpy/scipy/sklearn/dask (~150-250 MB each),
#   plus the client's own copy of the expression subset:
#     unpaired: 20,059 cells x 2,852 genes float32 = 0.23 GB
#   plus the scheduler and the final edge DataFrame (659,314 edges unpaired,
#   ~751k for an RNA-only baseline) held on the client at the end.
#
#   Rule of thumb:            --mem  ~=  (workers x 0.30 GB)
#                                     +  0.5 GB x (n_cells x n_cols / 1e9)  [client subset]
#                                     +  the edge-frame and scheduler overhead
#   =>  8 workers -> 32G is fine (32G is ample)
#      16 workers -> 64G  (this default; comfortable)
#      20 workers -> 96G  (max per PACE node)
#
#   If you hit OOM, in order: (1) lower GRN_WORKERS, (2) lower
#   GRN_THREADS_PER_WORKER.  Do NOT raise --mem past 96G for 20 workers --
#   memory scales with worker COUNT because of the broadcast, not with CPUs.
#
# TIME: the published unpaired run (20,059 x 2,852, 4 workers) took 2 h 10 m.
# Scaling measured on the real matrix (2,000 cells x 400 genes, seed 666,
# identical 21,998 edges at every setting):
#     1 worker 92.0s | 2 workers 47.2s | 4 workers 25.2s | 8 workers 18.6s
# Work is split per TARGET GENE, so the last gene to finish bounds the wall
# time and gains flatten past ~8-16 workers.  4 h covers 16 workers with margin.
# ---------------------------------------------------------------------------
#
# MULTI-NODE: set GRN_JOBQUEUE=1 to let this job submit its own worker jobs to
# SLURM (dask_jobqueue.SLURMCluster):
#     GRN_JOBQUEUE=1 GRN_NODE_WORKERS=2 GRN_SLURM_CORES=16 \
#       sbatch pace/slurm/04_grn.slurm
# That requires dask-jobqueue in the grn39 env (installed by 02_setup_envs.sh)
# and the submitting node being allowed to call sbatch.
#
#   sbatch --account=$PACE_ACCOUNT pace/slurm/04_grn.slurm
#   GRN_WORKERS=10 sbatch --account=$PACE_ACCOUNT --cpus-per-task=10 pace/slurm/04_grn.slurm
#   GRN_THREADS_PER_WORKER=2 GRN_WORKERS=8 ...   # 8 workers x 2 threads = 16 cores

set -euo pipefail
ENV_NAME="${GRN_ENV:-grn39}"          # Arboreto needs the python 3.9 env
source "$(dirname "$0")/_common.sh"

cd "$PACE_ROOT"
EXP="${PACE_EXPERIMENT:-unpaired}"

# One worker process per core by default; keep this == --cpus-per-task.
export GRN_WORKERS="${GRN_WORKERS:-$SLURM_CPUS_PER_TASK}"
export GRN_CPUS="${GRN_CPUS:-$SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$GRN_THREADS_PER_WORKER}"

echo "### GRN: experiment=$EXP workers=$GRN_WORKERS " \
     "threads_per_worker=${GRN_THREADS_PER_WORKER:-1} " \
     "mem_requested=$((SLURM_MEM_PER_NODE / 1024))G jobqueue=${GRN_JOBQUEUE:-0}"
echo "### (16 workers x 816 TFs x 20,059 cells x 4B = ~1.0 GB broadcast total; " \
     "client subset 0.23 GB)"

python scripts/06_run_grn.py --experiment "$EXP"

O="$PACE_ROOT/results/${EXP}_grn/grn_tfonly_reg"
echo
echo "=== run info ==="
cat "$O/run_info.txt"
echo
echo "=== outputs ==="
ls -la "$O"
