#!/bin/bash
# Preflight check -- run this after 02_setup_envs.sh / 03_patch_scsaga.py /
# 01_download_data.sh and BEFORE submitting, to catch configuration mistakes
# without burning queue time.
#
#   bash scripts/99_preflight.sh
#
# Exit code 0 = ready to submit.

set -uo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BUNDLE_DIR"
[ -f ./pace.env ] && source ./pace.env

PACE_ROOT="${PACE_ROOT:-$BUNDLE_DIR/workspace}"
PACE_DATA="${PACE_DATA:-$PACE_ROOT/data}"
PACE_RAW="${PACE_RAW:-$PACE_ROOT/raw}"
SCAGA_REPO="${SCAGA_REPO:-$PACE_ROOT/tools/scSAGA}"
EXP="${PACE_EXPERIMENT:-unpaired}"

FAIL=0
ok()   { printf '  \033[32mOK\033[0m   %s\n' "$*"; }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$*"; FAIL=1; }
warn() { printf '  \033[33mWARN\033[0m %s\n' "$*"; }
hdr()  { echo; echo "=== $*"; }

hdr "configuration"
echo "  bundle      : $BUNDLE_DIR"
echo "  PACE_ROOT   : $PACE_ROOT"
echo "  PACE_DATA   : $PACE_DATA"
echo "  PACE_RAW    : $PACE_RAW"
echo "  SCAGA_REPO  : $SCAGA_REPO"
echo "  experiment  : $EXP"
echo "  account     : ${PACE_ACCOUNT:-<unset>}   qos: ${PACE_QOS:-<unset>}"
echo "  GRN workers : ${GRN_WORKERS:-auto} x threads ${GRN_THREADS_PER_WORKER:-1}"

hdr "slurm"
if command -v sbatch >/dev/null 2>&1; then
  ok "sbatch present ($(command -v sbatch))"
else
  warn "sbatch not found -- expected off-cluster; the .slurm scripts cannot be tested here"
fi
[ -n "${PACE_ACCOUNT:-}" ] && ok "charge account set" \
  || warn "PACE_ACCOUNT unset -- sbatch will need --account, or jobs will be rejected"

hdr "conda / environments"
if command -v conda >/dev/null 2>&1; then
  ok "conda present"
  for e in "${SC_ENV:-scmint}" "${GRN_ENV:-grn39}"; do
    if conda env list | awk '{print $1}' | grep -qx "$e"; then
      PYV=$(conda run -n "$e" python -c 'import sys;print(".".join(map(str,sys.version_info[:3])))' 2>/dev/null)
      ok "env '$e' exists (python ${PYV:-?})"
      if [ "$e" = "${GRN_ENV:-grn39}" ]; then
        conda run -n "$e" python - <<'PY' >/dev/null 2>&1 \
          && ok "grn env: arboreto importable, click/_unicodefun present" \
          || bad "grn env: arboreto or click._unicodefun broken (need click==8.0.4)"
import arboreto, click
from click import _unicodefun
from arboreto.algo import grnboost2
PY
        conda run -n "$e" python -c "from dask_jobqueue import SLURMCluster" >/dev/null 2>&1 \
          && ok "grn env: dask-jobqueue present (multi-node GRN available)" \
          || warn "grn env: dask-jobqueue missing -> single-node GRN only"
      else
        conda run -n "$e" python - <<'PY' >/dev/null 2>&1 \
          && ok "scmint env: h5py/scipy/sklearn/ot/geosketch/torch importable" \
          || bad "scmint env: missing modules -- re-run envs/02_setup_envs.sh"
import h5py, scipy, sklearn, numpy, ot, geosketch, torch, yaml
PY
      fi
    else
      bad "env '$e' missing -- run envs/02_setup_envs.sh"
    fi
  done
else
  bad "conda not on PATH -- 'module load anaconda3' first"
fi

hdr "scSAGA checkout + SAVE-H patch"
if [ -d "$SCAGA_REPO/scmint" ]; then
  if grep -q "joint_embedding_H.npy" "$SCAGA_REPO/scmint/scsaga.py" 2>/dev/null; then
    ok "SAVE-H patch present in $SCAGA_REPO/scmint/scsaga.py"
  else
    bad "SAVE-H patch MISSING -- run: python envs/03_patch_scsaga.py"
  fi
  [ -f "$SCAGA_REPO/scmint/scsaga_saveH.py" ] \
    && ok "compatibility module scsaga_saveH.py present" \
    || warn "scmint/scsaga_saveH.py absent (only matters if something imports it)"
else
  bad "no scSAGA checkout at $SCAGA_REPO -- run: python envs/03_patch_scsaga.py"
fi

hdr "raw inputs (downloaded)"
# name:bytes pairs, checked with plain shell (no associative arrays, so this
# still runs under bash 3.2 as shipped on macOS)
check_raw() {
  f="$1"; want="$2"; p="$PACE_RAW/$f"
  if [ -f "$p" ]; then
    sz=$(wc -c < "$p" | tr -d ' ')
    if [ "$sz" = "$want" ]; then ok "$f ($sz bytes)"
    else warn "$f size $sz != expected $want"; fi
  else
    bad "$f missing -- run: bash scripts/01_download_data.sh"
  fi
}
check_raw "3k_multiome.h5" "38844318"
check_raw "10k_multiome.h5" "192125528"
check_raw "atac10k_v1.1.h5" "67041999"
check_raw "pbmc4k_filtered_gene_bc_matrices.tar.gz" "18423814"

hdr "prepared data (step 1)"
if [ "$EXP" = "unpaired" ]; then
  WANT="10k_rna 10k_atac atac10k_ext"
else
  WANT="10k_rna 10k_atac"
fi
for d in $WANT; do
  miss=""
  for f in counts.mtx barcodes.txt features.txt pca_50.txt; do
    [ -f "$PACE_DATA/$d/$f" ] || miss="$miss $f"
  done
  if [ -z "$miss" ]; then
    c=$(wc -l < "$PACE_DATA/$d/barcodes.txt")
    g=$(wc -l < "$PACE_DATA/$d/features.txt")
    ok "$d ($g features x $c cells)"
  else
    bad "$d incomplete (missing:$miss) -- run pace/slurm/01_preprocess.slurm"
  fi
done

hdr "repo reference / ground-truth files"
# The regulator/target lists and ground truth ship inside this repo.  Look in the
# obvious places rather than assuming a fixed relative depth.
REPO_DATA=""
for cand in "$BUNDLE_DIR/../data" "$BUNDLE_DIR/../../data" "$PACE_ROOT/repo/data" "$PACE_DATA"; do
  if [ -f "$cand/trrust_tf.txt" ] && [ -f "$cand/tf_only.txt" ]; then
    REPO_DATA="$(cd "$cand" && pwd)"; break
  fi
done
if [ -n "$REPO_DATA" ]; then
  ok "gene lists found in $REPO_DATA"
  echo "       (06_run_grn.py defaults to <repo>/data/trrust_tf.txt; override with TRRUST_FILE / TF_ONLY_FILE)"
  for f in ground_truth/PBMC-TRRUST.csv ground_truth/PBMC-Blood.csv; do
    [ -f "$REPO_DATA/$f" ] && ok "  $f" || bad "  $f missing (needed by 07_evaluate_grn.py)"
  done
else
  bad "trrust_tf.txt / tf_only.txt not found near the repo checkout"
  echo "       these ship in pbmc-multiome-grn-inference/data/ -- make sure you cloned the repo"
fi

hdr "storage"
if [ -d "$PACE_ROOT" ]; then
  avail=$(df -Pk "$PACE_ROOT" | awk 'NR==2{printf "%.1f", $4/1048576}')
  echo "  $PACE_ROOT : ${avail} GB free"
  awk -v a="$avail" 'BEGIN{exit !(a<40)}' \
    && warn "less than ~40 GB free; the pipeline needs roughly 3 GB for the all-cells matrix + working copies" \
    || ok "enough space"
else
  warn "$PACE_ROOT does not exist yet (will be created on submit)"
fi
case "$PACE_ROOT" in
  */home/*|$HOME) warn "PACE_ROOT is under home -- 20 GB quota; prefer ~/scratch" ;;
esac

hdr "result"
if [ "$FAIL" = "0" ]; then
  printf '  \033[32mREADY\033[0m -- submit with:  bash pace/submit_all.sh %s\n' "$EXP"
  exit 0
else
  printf '  \033[31mNOT READY\033[0m -- fix the FAIL lines above\n'
  exit 1
fi
