#!/bin/bash
# Create the three conda environments the pipeline needs.  Run once.
#
#   bash envs/setup_envs.sh
#
#   scmint   py3.12  integration (scSAGA), reverse-imputeKNN, evaluation
#   scement  py3.11  SCEMENT combined RNA reference (needs anndata/scanpy)
#   grn39    py3.9   Arboreto GRNBoost2
#
# WHY THREE: the three steps have mutually incompatible requirements.
#   * arboreto 0.1.6 is incompatible with modern numpy/dask -- it needs the
#     2021 numpy 1.21 / dask 2021.10 stack.
#   * distributed 2021.10 imports click._unicodefun, which click 8.1 removed:
#     with click >= 8.1 the dask-worker/dask-scheduler CLIs crash with
#       ImportError: cannot import name '_unicodefun' from 'click'
#     so grn39 pins click==8.0.4 (verified to ship _unicodefun).
#   * SCEMENT needs anndata/scanpy, kept out of the main env on purpose.
set -euo pipefail

SC_ENV="${SC_ENV:-scmint}"
SCE_ENV="${SCE_ENV:-scement}"
GRN_ENV="${GRN_ENV:-grn39}"
FORCE="${FORCE:-0}"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found -- on PACE run:  module load anaconda3" >&2
  module load anaconda3 2>/dev/null || true
fi
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda unavailable" >&2; exit 1; }

have_env() { conda env list | awk '{print $1}' | grep -qx "$1"; }
drop_env() { [ "$FORCE" = "1" ] && conda env remove -y -n "$1" 2>/dev/null || true; }

echo "=== env $SC_ENV (py3.12: scSAGA, imputation, evaluation)"
if have_env "$SC_ENV" && [ "$FORCE" != "1" ]; then echo "    exists, skipping"
else
  drop_env "$SC_ENV"
  conda create -y -n "$SC_ENV" -c conda-forge python=3.12 pip
fi
# CPU-only torch: scSAGA runs on CPU, this avoids pulling the CUDA wheels
conda run -n "$SC_ENV" pip -q install --index-url https://download.pytorch.org/whl/cpu torch
conda run -n "$SC_ENV" pip -q install \
  "numpy>=2,<3" "scipy>=1.13" "scikit-learn>=1.4" "pandas>=2.2" \
  "h5py>=3.10" pyyaml matplotlib joblib "pot>=0.9.4" geosketch

echo "=== env $SCE_ENV (py3.11: SCEMENT reference)"
if have_env "$SCE_ENV" && [ "$FORCE" != "1" ]; then echo "    exists, skipping"
else
  drop_env "$SCE_ENV"
  conda create -y -n "$SCE_ENV" -c conda-forge python=3.11 pip
fi
conda run -n "$SCE_ENV" pip -q install \
  "numpy>=2,<3" "scipy>=1.14" "pandas>=2.2" scikit-learn \
  anndata scanpy psutil joblib pyyaml

echo "=== env $GRN_ENV (py3.9: Arboreto GRNBoost2)"
if have_env "$GRN_ENV" && [ "$FORCE" != "1" ]; then echo "    exists, skipping"
else
  drop_env "$GRN_ENV"
  conda create -y -n "$GRN_ENV" -c conda-forge \
    "python=3.9" "numpy=1.21" "pandas=1.4" "scipy=1.9" "scikit-learn=1.1" \
    "dask=2021.10" "distributed=2021.10" "click=8.0.4" joblib pip
fi
conda run -n "$GRN_ENV" pip -q install arboreto==0.1.6 "dask-jobqueue==0.7.2"

echo
echo "=== verifying"
conda run -n "$SC_ENV" python - <<'PY'
import numpy, scipy, sklearn, h5py, ot, torch
print(f'  {SC_ENV}: numpy {numpy.__version__} scipy {scipy.__version__} '
      f'sklearn {sklearn.__version__} h5py {h5py.__version__} torch {torch.__version__}')
PY
conda run -n "$SCE_ENV" python - <<'PY'
import anndata, scanpy, numpy
print(f'  scement: anndata {anndata.__version__} scanpy {scanpy.__version__} numpy {numpy.__version__}')
PY
conda run -n "$GRN_ENV" python - <<'PY'
import numpy, dask, distributed, click
from arboreto.algo import grnboost2
print(f'  grn39: numpy {numpy.__version__} dask {dask.__version__} '
      f'distributed {distributed.__version__} click {click.__version__}')
try:
    from click import _unicodefun
    print('  grn39: click._unicodefun OK (dask CLI works)')
except ImportError:
    print('  grn39: WARN click._unicodefun missing -- pin click==8.0.4')
print('  grn39: arboreto importable')
PY

echo
echo "done.  use:  conda activate $SC_ENV | $SCE_ENV | $GRN_ENV"
