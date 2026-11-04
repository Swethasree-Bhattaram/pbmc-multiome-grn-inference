#!/bin/bash
# =============================================================================
#  Environment for scripts/stepwise_eval.py  (MultiVI + scSAGA, one env)
# =============================================================================
#   bash setup_stepwise.sh                # conda env "scgrn", adds muon+scvi-tools
#   ENV_NAME=mygrn bash setup_stepwise.sh
#   USE_VENV=1 bash setup_stepwise.sh     # plain ./.venv instead of conda
#   SKIP_BASE=1 bash setup_stepwise.sh    # only the MultiVI part (run.py env exists)
#   SMOKE=1 bash setup_stepwise.sh        # also self-test the metrics + build smoke data
#
#  WHAT IT ADDS
#  run.py's env (setup.sh) already has numpy/scipy/sklearn/pandas/scanpy/anndata,
#  arboreto + dask, and scSAGA (which pulls torch and scib).  The step-wise
#  experiment needs two more packages, and only for the MultiVI arm:
#
#      muon          the MuData container MultiVI is set up on
#      scvi-tools    MULTIVI itself
#
#  Everything else -- scSAGA, the reverse-imputeKNN, GRNBoost2, the metrics -- runs
#  on what setup.sh installed.  Versions this was verified against: scvi-tools
#  1.4.2, muon 0.1.9, torch 2.15, numpy 2.4, sklearn 1.9.
#
#  ON PACE: run this on a LOGIN NODE (compute nodes have no network).  Home is
#  quota'd at 10 GB, and three things default there even when the env lives on
#  scratch -- TMPDIR, pip's cache and conda's package dir.  Point all three at
#  scratch first, or torch alone will hit Errno 122:
#
#      export TMPDIR=/storage/scratch1/<n>/<user>/tmp
#      export PIP_CACHE_DIR=$TMPDIR/pip
#      export CONDA_PKGS_DIRS=$TMPDIR/pkgs
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
ENV_NAME="${ENV_NAME:-scgrn}"
USE_VENV="${USE_VENV:-0}"
_PY_GIVEN=""
SKIP_BASE="${SKIP_BASE:-0}"
SMOKE="${SMOKE:-0}"

# --- 1. the base env (integration + GRN + metrics) ---------------------------
# A pre-set PY wins: lets you point at an existing interpreter (a micromamba env,
# a venv outside the repo) without letting this script create one.
if [ -n "${PY:-}" ]; then
  _PY_GIVEN=1
  echo "=== using PY from the environment: $PY"
  "$PY" -c 'import sys; assert sys.version_info >= (3,10), sys.version' \
    || { echo "ERROR: PY must be python >= 3.10" >&2; exit 1; }
elif [ "$SKIP_BASE" = "1" ]; then
  echo "=== skipping setup.sh (SKIP_BASE=1)"
  if [ "$USE_VENV" = "1" ]; then PY="$ROOT/.venv/bin/python"
  else
    command -v conda >/dev/null 2>&1 || module load anaconda3 2>/dev/null || true
    source "$(conda info --base)/etc/profile.d/conda.sh"
    PY="$(conda run -n "$ENV_NAME" python -c 'import sys;print(sys.executable)')"
  fi
else
  USE_VENV="$USE_VENV" ENV_NAME="$ENV_NAME" bash "$ROOT/setup.sh"
  if [ "$USE_VENV" = "1" ]; then PY="$ROOT/.venv/bin/python"
  else
    command -v conda >/dev/null 2>&1 || module load anaconda3 2>/dev/null || true
    source "$(conda info --base)/etc/profile.d/conda.sh"
    PY="$(conda run -n "$ENV_NAME" python -c 'import sys;print(sys.executable)')"
  fi
fi
echo "    python: $PY"

# --- 2. MultiVI -----------------------------------------------------------------
echo "=== installing muon + scvi-tools (MultiVI arm)"
# a `uv venv` has no pip; prefer pip, then uv --python, then bootstrap pip
if "$PY" -m pip --version >/dev/null 2>&1; then
  PIP_INSTALL=("$PY" -m pip install -q)
elif command -v uv >/dev/null 2>&1; then
  PIP_INSTALL=(uv pip install -q --python "$PY")
elif "$PY" -m ensurepip --upgrade >/dev/null 2>&1; then
  PIP_INSTALL=("$PY" -m pip install -q)
else
  echo "ERROR: no pip and no uv for $PY" >&2; exit 1
fi
"${PIP_INSTALL[@]}" muon scvi-tools

# --- 3. verify -------------------------------------------------------------------
echo "=== verifying"
"$PY" - <<'PY'
import sys, importlib
mods = ["numpy", "scipy", "sklearn", "pandas", "anndata", "scanpy",
        "torch", "muon", "scvi", "arboreto", "dask", "distributed"]
bad = []
for m in mods:
    try:
        mod = importlib.import_module(m)
        print(f"  {m:12s} {getattr(mod, '__version__', '?')}")
    except Exception as exc:                                   # noqa: BLE001
        bad.append(m)
        print(f"  {m:12s} MISSING ({type(exc).__name__}: {exc})")
try:
    from scmint.scsaga import Saga
    print("  scSAGA       importable")
except Exception as exc:                                       # noqa: BLE001
    bad.append("scsaga")
    print(f"  scSAGA       MISSING ({exc})")
try:
    from scvi.model import MULTIVI
    import inspect
    assert hasattr(MULTIVI, "setup_mudata"), "MULTIVI.setup_mudata missing"
    print(f"  MULTIVI      importable  {inspect.signature(MULTIVI.__init__)}"
          .replace("(self, adata", "(")[:110])
except Exception as exc:                                       # noqa: BLE001
    bad.append("MULTIVI")
    print(f"  MULTIVI      MISSING ({exc})")
try:
    from arboreto.algo import grnboost2                    # noqa: F401
    import grn_compat
    grn_compat.apply()
    print("  arboreto     importable (+ modern-dask patch applied)")
except Exception as exc:                                       # noqa: BLE001
    bad.append("arboreto")
    print(f"  arboreto     FAILED ({exc})")
print("  python      ", sys.version.split()[0])
if bad:
    raise SystemExit(f"FAILED to import: {bad}")
print("  all imports OK")
PY

# --- 4. optional: prove the metrics and the plumbing on synthetic data -----------
if [ "$SMOKE" = "1" ]; then
  echo "=== metrics self-test (synthetic data with known answers)"
  "$PY" scripts/stepwise_metrics.py --self-test
  echo "=== building a tiny paired dataset under results/smoke"
  "$PY" scripts/make_smoke_data.py --out results/smoke/data
  echo "    then:  $PY scripts/stepwise_eval.py \\"
  echo "             --config results/smoke/config_smoke.yml --epochs 100 \\"
  echo "             --seeds 666 --workers 4 --out results/smoke/results"
fi

echo
if [ "$USE_VENV" = "1" ] || [ -n "${_PY_GIVEN:-}" ]; then
  echo "done.  run:  $PY scripts/stepwise_eval.py --check"
else
  echo "done.  run:  conda activate $ENV_NAME && python scripts/stepwise_eval.py --check"
fi
