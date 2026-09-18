#!/bin/bash
#SBATCH --job-name=grn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#
# =============================================================================
#  RUN THE PIPELINE ON PACE PHOENIX
# =============================================================================
#  Nothing in this file is tied to a particular path: it runs whatever
#  config.yml is next to it, wherever you cloned the repo.  Submit from the
#  repo directory, or point it there:
#
#      cd /storage/.../my-repo
#      sbatch --account=<acct> --qos=inferno run_grn.slurm
#      sbatch --account=<acct> --qos=inferno --chdir=/path/to/repo run_grn.slurm
#
#  Add --chdir / --mem / --time on the sbatch line to override the defaults
#  above.  --account and --qos are intentionally NOT baked in, so you can use
#  whichever allocation you like.
#
#  ONE-TIME SETUP (on a login node -- compute nodes usually have no network):
#      cd <repo> && bash setup.sh
#  Then this job just runs.
#
#  USEFUL OVERRIDES
#      EXPERIMENTS="ab" sbatch ... run_grn.slurm     # only one experiment
#      ONLY=grn,evaluate sbatch ... run_grn.slurm    # reuse earlier steps
#      GRN_WORKERS=8 sbatch --cpus-per-task=8 --mem=32G ... run_grn.slurm
# =============================================================================
set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$REPO"
CONFIG="${CONFIG:-$REPO/config.yml}"
[ -f "$CONFIG" ] || { echo "ERROR: no config at $CONFIG" >&2; exit 1; }

# --- environment ------------------------------------------------------------
# setup.sh created this once.  Adjust ENV_NAME if you named it differently.
ENV_NAME="${ENV_NAME:-scgrn}"
if ! command -v conda >/dev/null 2>&1; then
  module load anaconda3 2>/dev/null || true
fi
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
  PY="$(python -c 'import sys;print(sys.executable)')"
else
  # fall back to a venv built with USE_VENV=1
  PY="$REPO/.venv/bin/python"
  [ -x "$PY" ] || { echo "ERROR: neither conda env '$ENV_NAME' nor $PY" >&2
                    echo "       run: bash setup.sh   (or USE_VENV=1 bash setup.sh)" >&2
                    exit 1; }
fi

# --- GRN sizing -------------------------------------------------------------
# GRNBoost2 broadcasts the regulator matrix to every dask worker, so memory
# scales with WORKER COUNT, not with CPUs:  ~0.35 GB per worker + ~2 GB.
#    8 workers -> 32G      16 workers -> 64G      20 workers -> 96G
export GRN_WORKERS="${GRN_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}"
export GRN_THREADS_PER_WORKER="${GRN_THREADS_PER_WORKER:-1}"
export OMP_NUM_THREADS="$GRN_THREADS_PER_WORKER"
EXPERIMENTS="${EXPERIMENTS:-}"
ONLY="${ONLY:-integrate,impute,grn,evaluate}"

echo "==================== GRN PIPELINE ===================="
echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-<none>}  account=${SLURM_JOB_ACCOUNT:-none} qos=${SLURM_JOB_QOS:-none}"
echo "started     : $(date -Is)"
echo "repo        : $REPO"
echo "config      : $CONFIG"
echo "python      : $PY"
echo "cpus / mem  : ${SLURM_CPUS_PER_TASK:-?} / ${SLURM_MEM_PER_NODE:-?}MB"
echo "workers     : $GRN_WORKERS x $GRN_THREADS_PER_WORKER threads"
echo "steps       : $ONLY"
echo "experiments : ${EXPERIMENTS:-<all in config>}"
echo "======================================================"

ARGS=(--config "$CONFIG" --only "$ONLY")
for e in $EXPERIMENTS; do ARGS+=(-e "$e"); done

"$PY" run.py "${ARGS[@]}"

echo
echo "==================== RESULTS ========================="
find "$REPO/results" -name network.tsv -o -name evaluation_summary.txt 2>/dev/null | sort | sed 's/^/  /'
echo "finished    : $(date -Is)"
