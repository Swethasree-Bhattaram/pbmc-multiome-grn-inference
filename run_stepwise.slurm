#!/bin/bash
#SBATCH --job-name=stepwise
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#
# =============================================================================
#  STEP-WISE EXPERIMENT ON PACE PHOENIX: MultiVI vs scSAGA, every step measured
# =============================================================================
#  Runs scripts/stepwise_eval.py on the config.yml next to this script.  Nothing
#  is tied to a path: submit from the repo, or point at it.
#
#      cd /storage/.../pbmc-multiome-grn-inference
#      sbatch --account=<acct> --qos=inferno run_stepwise.slurm
#      sbatch --account=<acct> --qos=inferno --chdir=/path/to/repo run_stepwise.slurm
#
#  ONE-TIME SETUP on a login node (compute nodes have no network):
#      bash setup_stepwise.sh
#      python scripts/stepwise_eval.py --check
#
#  WHAT IT RUNS (in this order, per method, per experiment)
#      integrate  scSAGA joint embedding / MultiVI fit -> H (+ model)
#      impute     reverse-imputeKNN / MultiVI decode -> all-cells matrix + mask
#                 plus the frozen held-out arm (a second MultiVI fit for MultiVI)
#      grn        GRNBoost2, one network per seed
#      metrics    every metric in workflow.pptx -> metrics.json + report.md
#
#  OVERRIDES (all optional)
#      EXPERIMENTS="paired" sbatch ... run_stepwise.slurm     # one experiment
#      METHODS="scsaga"     sbatch ... run_stepwise.slurm     # one workflow
#      STEPS="metrics"      sbatch ... run_stepwise.slurm     # reuse earlier steps
#      EPOCHS=100           sbatch ... run_stepwise.slurm     # MultiVI max epochs
#      SEEDS="666"          sbatch ... run_stepwise.slurm     # one network, no stability
#      GRN_WORKERS=8 sbatch --cpus-per-task=8 --mem=32G ... run_stepwise.slurm
#      FORCE=1              sbatch ... run_stepwise.slurm     # redo existing networks
#      MULTIVI_THREADS=16   sbatch ... run_stepwise.slurm     # torch threads for MultiVI
#
#  THREADS: OMP_NUM_THREADS is pinned to GRN_THREADS_PER_WORKER because GRNBoost2
#  runs one thread per dask worker.  That is wrong for a CPU-only MultiVI fit, so
#  stepwise_eval sets torch's intra-op pool itself, to all cores minus one unless
#  MULTIVI_THREADS says otherwise (measured at 36,601 genes + 89,796 peaks:
#  7.1 s/epoch single-threaded vs 3.2 s/epoch with 8 threads).
#
#  MEMORY: GRNBoost2 broadcasts the regulator matrix to every dask worker, so
#  memory scales with WORKER COUNT, not with CPUs (8 -> 32G, 16 -> 64G, 20 -> 96G),
#  and MultiVI holds the RNA + ATAC count matrices in memory on top of that.  64G
#  at 16 workers is the default here; raise --mem if you also run the biggest
#  experiment with several hundred thousand peaks.
#
#  TIME: the 3-seed GRN loop is the expensive part and it runs twice (once per
#  method).  For a first pass use SEEDS="666" to get one network per method, then
#  re-submit with STEPS=grn SEEDS="667,668" once you know the input is right.
# =============================================================================
set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$PWD}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -f "$REPO/config.yml" ] && [ -f "$SCRIPT_DIR/config.yml" ]; then
  REPO="$SCRIPT_DIR"
fi
cd "$REPO"
CONFIG="${CONFIG:-$REPO/config.yml}"
[ -f "$CONFIG" ] || { echo "ERROR: no config at $CONFIG" >&2; exit 1; }

# --- environment ------------------------------------------------------------
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
  PY="$REPO/.venv/bin/python"
  [ -x "$PY" ] || { echo "ERROR: neither conda env '$ENV_NAME' nor $PY" >&2
                    echo "       run: bash setup_stepwise.sh" >&2
                    exit 1; }
fi

# --- knobs ------------------------------------------------------------------
export GRN_WORKERS="${GRN_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}"
export GRN_THREADS_PER_WORKER="${GRN_THREADS_PER_WORKER:-1}"
export GRN_SEEDS="${SEEDS:-666,667,668}"
export OMP_NUM_THREADS="$GRN_THREADS_PER_WORKER"
export MKL_NUM_THREADS="$GRN_THREADS_PER_WORKER"
METHODS="${METHODS:-multivi,scsaga}"
STEPS="${STEPS:-integrate,impute,grn,metrics}"
EXPERIMENTS="${EXPERIMENTS:-}"
EPOCHS="${EPOCHS:-200}"
ARGS=(--config "$CONFIG" --only "$STEPS" --methods "$METHODS" --epochs "$EPOCHS")
for e in $EXPERIMENTS; do ARGS+=(-e "$e"); done
[ "${FORCE:-0}" = "1" ] && ARGS+=(--force)

echo "==================== STEP-WISE EXPERIMENT ===================="
echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-<none>}  account=${SLURM_JOB_ACCOUNT:-none} qos=${SLURM_JOB_QOS:-none}"
echo "started     : $(date '+%Y-%m-%dT%H:%M:%S%z')"
echo "repo        : $REPO"
echo "config      : $CONFIG"
echo "python      : $PY"
echo "cpus / mem  : ${SLURM_CPUS_PER_TASK:-?} / ${SLURM_MEM_PER_NODE:-?}MB"
echo "methods     : $METHODS"
echo "steps       : $STEPS"
echo "experiments : ${EXPERIMENTS:-<all in config>}"
echo "epochs/seeds: $EPOCHS / $GRN_SEEDS"
echo "workers     : $GRN_WORKERS x $GRN_THREADS_PER_WORKER threads"
echo "=============================================================="

"$PY" scripts/stepwise_eval.py "${ARGS[@]}"

echo
echo "==================== RESULTS ================================="
find "$REPO/results/stepwise" -name network.tsv -o -name metrics.json -o -name report.md 2>/dev/null \
  | sort | sed 's/^/  /'
echo
for r in "$REPO"/results/stepwise/*/report.md; do
  [ -f "$r" ] || continue
  echo "----- $r"
  sed -n '1,60p' "$r"
done
echo "finished    : $(date '+%Y-%m-%dT%H:%M:%S%z')"
