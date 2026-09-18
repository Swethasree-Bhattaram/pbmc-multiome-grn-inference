#!/bin/bash
# Create the two conda environments the pipeline needs, on PACE Phoenix.
#
#   scmint   python 3.12  -- h5 splitting, scSAGA integration, reverse-imputeKNN,
#                            QC, evaluation, reports
#   grn39    python 3.9   -- Arboreto GRNBoost2 (needs the old dask/distributed stack)
#
# Run on a LOGIN node (installs only, no compute):
#   bash envs/02_setup_envs.sh
#
# Idempotent: existing envs are skipped unless FORCE=1.
#
# WHY THE PINS
#   * arboreto 0.1.6 is incompatible with modern numpy/dask; it needs
#     python 3.9 + numpy 1.21 + dask/distributed 2021.10.
#   * click MUST be < 8.1 in the py3.9 env.  distributed 2021.10 imports
#     `click._unicodefun`, which click 8.1 removed -- with click >= 8.1 the
#     `dask-worker` / `dask-scheduler` command line tools crash on import:
#       ImportError: cannot import name '_unicodefun' from 'click'
#     Pinning click==8.0.4 fixes it (verified: click 8.0.4 ships _unicodefun).
#   * Torch is installed CPU-only; the pipeline runs scSAGA on CPU.
set -euo pipefail

ROOT="${PACE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
FORCE="${FORCE:-0}"
SCAGA_REPO="${SCAGA_REPO:-$ROOT/tools/scSAGA}"
SC_ENV="${SC_ENV:-scmint}"
GRN_ENV="${GRN_ENV:-grn39}"

echo "PACE_ROOT   = $ROOT"
echo "scSAGA repo = $SCAGA_REPO"
echo "scSAGA env  = $SC_ENV"
echo "arboreto env= $GRN_ENV"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found. On PACE:  module load anaconda3" >&2
  echo "(loading it for this shell)" >&2
  module load anaconda3 2>/dev/null || true
fi
command -v conda >/dev/null 2>&1 || {
  echo "ERROR: conda still unavailable. Run 'module load anaconda3' then re-run." >&2
  exit 1
}

# ---------------------------------------------------------------- scSAGA env
if conda env list | awk '{print $1}' | grep -qx "$SC_ENV" && [ "$FORCE" != "1" ]; then
  echo "=== env '$SC_ENV' exists -- skipping (FORCE=1 to rebuild)"
else
  [ "$FORCE" = "1" ] && conda env remove -y -n "$SC_ENV" 2>/dev/null || true
  echo "=== creating env '$SC_ENV' (python 3.12)"
  conda create -y -n "$SC_ENV" -c conda-forge python=3.12 pip
fi
SC_PY="$(conda run -n "$SC_ENV" python -c 'import sys; print(sys.executable)')"
echo "    python: $SC_PY"
conda run -n "$SC_ENV" pip -q install --upgrade pip

# CPU-only torch first, so the CUDA wheels are never pulled in.
conda run -n "$SC_ENV" pip -q install \
  --index-url https://download.pytorch.org/whl/cpu torch
conda run -n "$SC_ENV" pip -q install \
  "numpy>=2,<3" "scipy>=1.13" "scikit-learn>=1.4" "pandas>=2.2" \
  "h5py>=3.10" pyyaml matplotlib joblib pydantic \
  "pot>=0.9.4" geosketch

# --------------------------------------------------------------- arboreto env
if conda env list | awk '{print $1}' | grep -qx "$GRN_ENV" && [ "$FORCE" != "1" ]; then
  echo "=== env '$GRN_ENV' exists -- skipping (FORCE=1 to rebuild)"
else
  [ "$FORCE" = "1" ] && conda env remove -y -n "$GRN_ENV" 2>/dev/null || true
  echo "=== creating env '$GRN_ENV' (python 3.9, old dask stack)"
  conda create -y -n "$GRN_ENV" -c conda-forge \
    "python=3.9" "numpy=1.21" "pandas=1.4" "scipy=1.9" "scikit-learn=1.1" \
    "dask=2021.10" "distributed=2021.10" "click=8.0.4" joblib pip
fi
GRN_PY="$(conda run -n "$GRN_ENV" python -c 'import sys; print(sys.executable)')"
echo "    python: $GRN_PY"
conda run -n "$GRN_ENV" pip -q install arboreto==0.1.6 "dask-jobqueue==0.7.2"

# ---------------------------------------------------------------- verify envs
echo
echo "=== verifying env '$SC_ENV'"
conda run -n "$SC_ENV" python - <<'PY'
import importlib
mods = ['numpy', 'scipy', 'sklearn', 'pandas', 'h5py', 'yaml', 'matplotlib',
        'joblib', 'pydantic', 'ot', 'geosketch', 'torch']
for m in mods:
    try:
        mod = importlib.import_module(m)
        print(f'  OK  {m:12s} {getattr(mod, "__version__", "")}')
    except Exception as e:
        print(f'  FAIL {m}: {e}')
PY

echo
echo "=== verifying env '$GRN_ENV'"
conda run -n "$GRN_ENV" python - <<'PY'
import numpy, pandas, scipy, sklearn, dask, distributed, click, arboreto
print(f'  numpy {numpy.__version__}  pandas {pandas.__version__}  scipy {scipy.__version__}')
print(f'  sklearn {sklearn.__version__}  dask {dask.__version__}  distributed {distributed.__version__}')
print(f'  click {click.__version__}')
try:
    from click import _unicodefun
    print('  OK  click._unicodefun present -> dask CLI will work')
except ImportError:
    print('  FAIL click._unicodefun missing -> pin click==8.0.4')
from arboreto.algo import grnboost2
print('  OK  arboreto.algo.grnboost2 importable')
try:
    from dask_jobqueue import SLURMCluster
    print('  OK  dask_jobqueue.SLURMCluster importable (multi-node GRN available)')
except ImportError as e:
    print(f'  WARN dask_jobqueue unavailable ({e}) -> single-node GRN only')
PY

echo
echo "Environments ready."
echo "  scSAGA/preprocessing : conda activate $SC_ENV"
echo "  Arboreto GRN         : conda activate $GRN_ENV"
