# pace-phoenix

Run bundle for the scSAGA → reverse-imputeKNN → Arboreto GRNBoost2 pipeline on
the PACE **Phoenix** cluster.

**Start here: [README_PACE.md](README_PACE.md)** — quick start, dataset table with
MD5s, SLURM resource justification, the required scSAGA patch, and measured
Arboreto scaling numbers.

```
pace.env                  all site settings (workspace, account, qos, workers)
requirements.txt          verified package versions per env
README_PACE.md            full guide
envs/                     conda env setup + scSAGA clone/patch
scripts/                  pipeline steps, portable, env-driven
pace/submit_all.sh        submit the whole chain (afterok dependencies)
pace/slurm/*.slurm        one script per step
```

TL;DR on a Phoenix login node:

```bash
cd pace-phoenix
export PACE_ROOT=$HOME/scratch/pbmc-grn PACE_ACCOUNT=<your-MAM-account>
source pace.env
module load anaconda3
bash envs/02_setup_envs.sh          # conda envs scmint + grn39
python envs/03_patch_scsaga.py      # REQUIRED scSAGA SAVE-H patch
bash scripts/01_download_data.sh    # 4 raw 10x datasets, MD5-verified
bash pace/submit_all.sh unpaired    # download -> pre -> integrate -> impute -> GRN -> eval
```
