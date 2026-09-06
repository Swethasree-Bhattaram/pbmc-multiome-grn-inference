#!/bin/bash
# Run the three Arboreto GRN experiments SEQUENTIALLY with TRRUST-only target
# genes (full-10k), under `caffeinate` so the Mac cannot sleep and kill dask.
# Reuses existing all-cells matrices; writes to results/<exp>_10k/grn_trrust/.
cd /Volumes/samsung_ssd/tmp/pbmc-full10k-grn
PY39=/Volumes/samsung_ssd/tmp/scSAGA/.venv39/bin/python
PY=/Volumes/samsung_ssd/tmp/scSAGA/.venv/bin/python
export PBSC4K_ROOT=$(pwd)

caffeinate -dimsu -w $$ &
CAFF=$!
echo "caffeinate pid=$CAFF started $(date +%H:%M:%S)"

run_one() {
  exp=$1
  mkdir -p "results/${exp}_10k/grn_trrust"
  echo "=== START $exp arboreto $(date +%H:%M:%S) ==="
  ( ulimit -n 8192; exec "$PY39" scripts/run_arboreto_trrust_10k.py "$exp" ) \
       > "results/${exp}_10k/grn_trrust/arboreto_run.log" 2>&1
  echo "=== END $exp arboreto exit=$? $(date +%H:%M:%S) ==="
  echo "=== START $exp eval $(date +%H:%M:%S) ==="
  "$PY" scripts/evaluate_grn_trrust_10k.py "$exp" > "results/${exp}_10k/grn_trrust/evaluation_run.log" 2>&1
  echo "=== END $exp eval exit=$? $(date +%H:%M:%S) ==="
}

run_one expA
run_one expB1
run_one expB2

kill $CAFF 2>/dev/null
echo "ALL DONE $(date +%H:%M:%S)"
