#!/bin/bash
# Arboreto GRNBoost2 for the UNPAIRED 10k run:
#   RNA10k multiome (11,898 real) + external 10k ATAC v1.1 (8,161 imputed)
#   targets  = trrust_tf.txt genes present (2,827)
#   regulators = tf_only.txt TFs present (816)
# Run under `caffeinate` so the Mac cannot sleep and kill the dask workers.
set -u
cd /Volumes/samsung_ssd/tmp/pbmc-10k-unpaired-grn
PY39=/Volumes/samsung_ssd/tmp/scSAGA/.venv39/bin/python
export UNPAIRED_ROOT=$(pwd)
export GRN_WORKERS=${GRN_WORKERS:-4}

caffeinate -dimsu -w $$ &
CAFF=$!
echo "caffeinate pid=$CAFF started $(date +%H:%M:%S)"

mkdir -p results/unpaired_grn/grn_tfonly_reg
echo "=== START arboreto unpaired (workers=$GRN_WORKERS) $(date +%H:%M:%S) ==="
( ulimit -n 8192; exec "$PY39" scripts/run_arboreto_unpaired.py ) \
     > results/unpaired_grn/grn_tfonly_reg/arboreto_run.log 2>&1
RC=$?
echo "=== END arboreto unpaired exit=$RC $(date +%H:%M:%S) ==="

kill $CAFF 2>/dev/null
echo "ALL DONE $(date +%H:%M:%S)"
