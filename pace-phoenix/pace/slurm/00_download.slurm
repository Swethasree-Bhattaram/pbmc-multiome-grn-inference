#!/bin/bash
#SBATCH --job-name=grn-download
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm_%x-%j.out
#SBATCH --error=logs/slurm_%x-%j.err
#
# STEP 0 -- download + verify the four raw 10x datasets (~316 MB total).
#
# You can also just run `bash scripts/01_download_data.sh` on a login node; that
# is usually faster than waiting in a queue. This script exists so the whole
# pipeline can be kicked off from one sbatch chain if you prefer.
#
#   sbatch pace/slurm/00_download.slurm
#   sbatch --account=GT-xxxx pace/slurm/00_download.slurm   # override the account

set -euo pipefail
source "$(dirname "$0")/_common.sh"

cd "$PACE_ROOT"
bash scripts/01_download_data.sh

echo
echo "raw inputs:"
ls -la "$PACE_RAW"
