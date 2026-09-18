#!/usr/bin/env python
"""Run scSAGA integration using the modified scsaga_saveH.py so the full 30-dim
joint embedding H is saved directly (no separate extract_embedding.py needed).

Uses the new_benchmark3-free uv venv at /Volumes/samsung_ssd/tmp/scSAGA/.venv
and the modified source at /Volumes/samsung_ssd/tmp/scSAGA/scmint/scsaga_saveH.py.
Writes results to results/scsaga_output.
"""
import sys, os, time
sys.path.insert(0, os.environ.get("SCAGA_REPO", os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "scSAGA"))))
from scmint.scsaga_saveH import main as saga_main

CFG = os.environ.get("SCSAGA_CONFIG", os.path.join(
    os.environ.get("PBMC3K_ROOT", os.path.expanduser("~/pbmc3k_analysis")),
    "scripts", "scsaga_config.yml"))
import yaml
with open(CFG) as f:
    raw_cfg = yaml.safe_load(f)

t0 = time.time()
saga_main(raw_cfg)
print(f"\n[RUNNER] scSAGA main() finished in {time.time()-t0:.1f}s")
