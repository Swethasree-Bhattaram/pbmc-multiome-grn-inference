#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/home/hcoda1/8/sbhattaram6/sc_multimodal_ot/slurm_grn-%j.out
#SBATCH --error=/storage/home/hcoda1/8/sbhattaram6/sc_multimodal_ot/slurm_grn-%j.err
#SBATCH --chdir=/storage/home/hcoda1/8/sbhattaram6/sc_multimodal_ot
#
# =============================================================================
#  RUN THE GRN PIPELINE -- PACE Phoenix
# =============================================================================
#  For every experiment in config.yml:
#     scSAGA integration  ->  reverse-imputeKNN (one per RNA reference)
#     ->  Arboreto GRNBoost2  ->  evaluation
#
#  Add --account / --qos / --mail if you want them baked in, e.g.
#     #SBATCH --account=<your-MAM-account>
#     #SBATCH --qos=inferno                 # not preempted; embers is
#     #SBATCH --mail-type=END,FAIL
#     #SBATCH --mail-user=<you>@gatech.edu
#  or pass them on the sbatch line instead (see below).
#
# -----------------------------------------------------------------------------
# RESOURCES -- why these numbers
#
#   --cpus-per-task=16   the GRN step runs 16 dask worker PROCESSES. scSAGA and
#                        the imputation are far lighter, so GRN sets the CPU.
#
#   --mem=64G            set by the GRN step. Arboreto BROADCASTS the TF matrix
#                        to every worker (arboreto/core.py,
#                        client.scatter(tf_matrix, broadcast=True)):
#
#                            816 TFs x n_cells x 4 bytes  per worker
#
#                        For 29,218 cells that is 95 MB/worker -> 1.5 GB across
#                        16 workers; for 20,059 cells it is 65 MB -> 1.0 GB.
#                        Add ~150-250 MB/worker for the interpreter stack.
#
#                        Sizing rule:   --mem  =  workers x 0.35 GB  +  ~2 GB
#                            8 workers -> 32G
#                           16 workers -> 64G   <-- this file
#                           20 workers -> 96G   (PACE node max is ~192 GB)
#
#                        Memory scales with worker COUNT, not with CPUs, because
#                        of that broadcast. The imputation step also peaks here:
#                        the all-cells matrix is 29,218 x 36,601 float32 = 4.0 GB.
#
#   --time=24:00:00      the GRN step dominates. Measured scaling on a real
#                        matrix (2,000 cells x 400 genes, seed 666, identical
#                        21,998 edges at every worker count):
#                            1 worker  92.0 s
#                            2 workers 47.2 s
#                            4 workers 25.2 s
#                            8 workers 18.6 s
#                        The published 4-worker unpaired run (20,059 x 2,852)
#                        took 2 h 10 m. At 16 workers expect roughly 40-60 min
#                        PER strategy; 5 strategies across all experiments is
#                        ~4-7 h. 24 h leaves comfortable headroom.
#
#                        Use --qos=inferno for anything this long -- embers is
#                        free but preemptible after a guaranteed first hour.
#
#   To shrink:  GRN_WORKERS=8  with  --cpus-per-task=8 --mem=32G  (about 2x time)
# -----------------------------------------------------------------------------
# USAGE
#   sbatch --account=<acct> --qos=inferno slurm_grn.slurm          # all experiments
#   EXPERIMENTS="ab" sbatch ... slurm_grn.slurm                    # just one
#   SKIP_PREP=1 sbatch ... slurm_grn.slurm                         # data already built
#   GRN_WORKERS=8 sbatch --cpus-per-task=8 --mem=32G ... slurm_grn.slurm
# =============================================================================

set -euo pipefail

# --- where everything lives --------------------------------------------------
# The repo root: config.yml, data/, results/.  Either cd here (--chdir above) or
# override.  On PACE put this on scratch, not home (home is 20 GB).
PACE_ROOT="${PACE_ROOT:-$PWD}"
CONFIG="${CONFIG:-$PACE_ROOT/config.yml}"

SC_ENV="${SC_ENV:-scmint}"
SCE_ENV="${SCE_ENV:-scement}"
GRN_ENV="${GRN_ENV:-grn39}"

# scSAGA checkout.  ANY clone works -- this pipeline calls the Saga class and
# saves the joint embedding itself, so no patched fork is needed:
#     git clone https://github.com/AluruLab/scSAGA.git tools/scSAGA
export SCAGA_REPO="${SCAGA_REPO:-$PACE_ROOT/tools/scSAGA}"

# GRN worker count -- keep equal to --cpus-per-task.
export GRN_WORKERS="${GRN_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}"
export GRN_THREADS_PER_WORKER="${GRN_THREADS_PER_WORKER:-1}"
export PIPELINE_ROOT="$PACE_ROOT"
export OMP_NUM_THREADS="${GRN_THREADS_PER_WORKER}"

# which experiments (blank = every experiment in config.yml)
EXPERIMENTS="${EXPERIMENTS:-}"

# skip data steps if you already have them
SKIP_PREP="${SKIP_PREP:-0}"
SKIP_DOWNLOAD="${SKIP_DOWNLOAD:-0}"

echo "==================== GRN PIPELINE ===================="
echo "host       : $(hostname)"
echo "job        : ${SLURM_JOB_ID:-<none>}  account=${SLURM_JOB_ACCOUNT:-none} qos=${SLURM_JOB_QOS:-none}"
echo "date       : $(date -Is)"
echo "PACE_ROOT  : $PACE_ROOT"
echo "config     : $CONFIG"
echo "cpus / mem : ${SLURM_CPUS_PER_TASK:-?} / $(( ${SLURM_MEM_PER_NODE:-0} / 1024 ))G"
echo "GRN        : $GRN_WORKERS workers x $GRN_THREADS_PER_WORKER threads"
echo "scSAGA     : $SCAGA_REPO"
echo "experiments: ${EXPERIMENTS:-<all in config>}"
echo "======================================================"

module load anaconda3
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
cd "$PACE_ROOT"
mkdir -p logs

# --- environments + scSAGA + the whole pipeline -------------------------------
# run_all.sh does: create envs -> install scSAGA -> validate -> integrate ->
# impute -> GRN -> evaluate (+ comparison when >1 RNA dataset).
SKIP_ENVS="${SKIP_ENVS:-1}"   # SLURM jobs shouldn't build envs; do that once first
EXPERIMENTS="$EXPERIMENTS" \
SKIP_ENVS="$SKIP_ENVS" \
GRN_WORKERS="$GRN_WORKERS" \
GRN_THREADS_PER_WORKER="$GRN_THREADS_PER_WORKER" \
PIPELINE_ROOT="$PACE_ROOT" \
  bash run_all.sh

echo; echo "### done $(date -Is)"
echo
echo "networks:"
find "$PACE_ROOT/results" -name network.tsv | sort | sed 's/^/  /'
echo "evaluations:"
find "$PACE_ROOT/results" -name evaluation_summary.txt | sort | sed 's/^/  /'
