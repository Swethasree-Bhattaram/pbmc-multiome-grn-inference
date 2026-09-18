#!/bin/bash
# =============================================================================
#  One environment for the whole pipeline.
# =============================================================================
#   bash setup.sh                    # conda env named scgrn (default)
#   ENV_NAME=mygrn bash setup.sh     # different name
#   USE_VENV=1 bash setup.sh         # plain venv at ./.venv, no conda
#
#  Run this once, on a login node on PACE (compute nodes usually have no
#  network), then submit run_grn.slurm.
#
#  WHY ONE ENV IS ENOUGH
#  scSAGA already requires python>=3.10 and depends on anndata/scanpy, so the
#  integration, the reference combination (ComBat) and the GRN step can share
#  one interpreter.  The only reason the old layout needed three envs was
#  arboreto's last release pinning ancient numpy/dask; run.py patches that one
#  incompatibility at import time instead.
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
ENV_NAME="${ENV_NAME:-scgrn}"
USE_VENV="${USE_VENV:-0}"

# --- 1. create the environment ----------------------------------------------
if [ "$USE_VENV" = "1" ]; then
  PY_BIN="${PY:-python3}"
  echo "=== venv at $ROOT/.venv  ($("$PY_BIN" -V 2>&1))"
  [ -d "$ROOT/.venv" ] || "$PY_BIN" -m venv "$ROOT/.venv"
  PY="$ROOT/.venv/bin/python"
  PIP="$ROOT/.venv/bin/pip"
else
  if ! command -v conda >/dev/null 2>&1; then
    echo "--- conda not on PATH, trying: module load anaconda3"
    module load anaconda3 2>/dev/null || true
  fi
  command -v conda >/dev/null 2>&1 || {
    echo "ERROR: no conda.  On PACE: module load anaconda3" >&2
    echo "       Or run with USE_VENV=1 to build a plain venv instead." >&2
    exit 1; }
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  echo "=== conda env $ENV_NAME"
  conda env list | awk '{print $1}' | grep -qx "$ENV_NAME" \
    || conda create -y -n "$ENV_NAME" -c conda-forge python=3.11 pip
  PY="$(conda run -n "$ENV_NAME" python -c 'import sys;print(sys.executable)')"
  PIP="$(conda run -n "$ENV_NAME" python -m pip --version >/dev/null && echo "$PY -m pip")"
fi

echo "    python: $PY"
"$PY" -c 'import sys; assert sys.version_info >= (3,10), sys.version; print("    version OK", sys.version.split()[0])'

# --- 2. dependencies --------------------------------------------------------
echo "=== installing requirements"
"$PY" -m pip install -q --upgrade pip
"$PY" -m pip install -q -r "$ROOT/requirements.txt"

# scSAGA needs torch; the CPU wheel avoids pulling multi-GB CUDA deps.
echo "=== installing torch (CPU)"
"$PY" -m pip install -q --index-url https://download.pytorch.org/whl/cpu torch

echo "=== installing scSAGA from GitHub"
"$PY" -m pip install -q "scsaga @ git+https://github.com/AluruLab/scSAGA.git"

# --- 3. verify --------------------------------------------------------------
echo "=== verifying"
"$PY" - <<'PY'
import sys
import numpy, scipy, sklearn, pandas, anndata, scanpy, torch
import arboreto, dask, distributed
from arboreto.algo import grnboost2
from scmint.scsaga import Saga
import grn_compat
grn_compat.apply()
print(f'  python      {sys.version.split()[0]}')
print(f'  numpy       {numpy.__version__}')
print(f'  scipy       {scipy.__version__}')
print(f'  sklearn     {sklearn.__version__}')
print(f'  pandas      {pandas.__version__}')
print(f'  torch       {torch.__version__}')
print(f'  anndata     {anndata.__version__}')
print(f'  scanpy      {scanpy.__version__}')
print(f'  dask        {dask.__version__}')
print(f'  distributed {distributed.__version__}')
print('  scSAGA      importable')
print('  arboreto    importable (+ modern-dask patch applied)')
PY

echo
if [ "$USE_VENV" = "1" ]; then
  echo "done.  run:  python run.py --config config.yml"
else
  echo "done.  run:  conda activate $ENV_NAME && python run.py --config config.yml"
  echo "  or from a script/SLURM:  $PY run.py --config config.yml"
fi
