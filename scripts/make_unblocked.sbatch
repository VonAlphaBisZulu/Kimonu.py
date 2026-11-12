#!/usr/bin/env bash
#SBATCH --job-name=mkunb
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=8:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts/unblocked/slurm_%j.log
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
export PYTHONPATH=$HOME/sd-scip-pricing
ls $BASE/scripts/sweep_unblocked/*_fixed.crn.json | xargs -P 16 -I{} bash -c '
  f={}; m=$(basename $f _fixed.crn.json)
  BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
  $HOME/miniconda3/envs/straindesign/bin/python $BASE/scripts/write_unblocked_crn.py \
    $f $BASE/scripts/unblocked/${m}_fixed.crn.json --solver scip \
    > $BASE/scripts/unblocked/${m}.log 2>&1 || echo "FAILED $m"
'
echo "ALL DONE $(date -Is)"
