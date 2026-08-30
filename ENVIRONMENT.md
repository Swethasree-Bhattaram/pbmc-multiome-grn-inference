# Environment setup

This project runs three distinct steps, each needing a specific Python stack.

## 1. Data preprocessing & scSAGA integration

Needs: `numpy`, `scipy`, `scikit-learn`, `h5py`, plus the AluruLab `scSAGA`
package (`scmint.scsaga`) and its deps (`pot`, `geosketch`, `torch`, `pydantic`).

The modified `scmint/scsaga_saveH.py` (which persists the full joint embedding H)
lives in the scSAGA checkout. The scripts import it and run on CPU.

```
pip install numpy scipy scikit-learn h5py scSAGA torch pot geosketch pydantic
```

## 2. SCEMENT combined-RNA reference (Experiment A)

SCEMENT (`AluruLab/scement`) is ComBat-style batch integration. The full C++
build is optional; the pure-Python `sct_sparse` path is used here. Requires
`numpy<2` (SCEMENT uses `np.float_`), `scipy`, `pandas`, `anndata`, `scanpy`,
`formulaic`, `patsy`, `psutil`.

A numpy2/scipy2-compatible copy of `sct_sparse` is vendored at
`scripts/multi-dataset-4x/scement_py/scement_sparse.py` (only `.A`->`.toarray()`
and `np.float_`->`np.float64` fixes applied). Use Python 3.11 with the above deps.

```
python3.11 -m venv .venv-scement
.venv-scement/bin/pip install "numpy<2" scipy pandas anndata scanpy formulaic patsy psutil
```

## 3. Arboreto GRNBoost2

Arboreto 0.1.6 is incompatible with modern numpy/dask. Use a dedicated old-stack
env (python 3.9, numpy 1.21, dask 2021.10, distributed 2021.10, arboreto 0.1.6).

```
conda create -y -n grn3 -c conda-forge "python=3.9" "numpy=1.21" "pandas=1.3" \
  "scipy=1.7" "scikit-learn=1.0" "dask=2021.10" "distributed=2021.10" joblib numba=0.55
~/.conda/envs/grn3/bin/pip install arboreto
```

**Critical:** on macOS every Arboreto script must be wrapped in
`if __name__ == '__main__':` (multiprocessing uses `spawn`).

## Setting the project root

All multi-dataset scripts resolve the repo root from the `PBSC4K_ROOT`
environment variable (fallback: two parents above the script). Set:

```
export PBSC4K_ROOT=/path/to/this/repo
```
