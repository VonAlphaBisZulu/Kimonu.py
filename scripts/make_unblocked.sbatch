#!/usr/bin/env bash
#SBATCH --job-name=mkunb
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=8:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_records/networks/unblocked/slurm_%j.log
# The networks of the paper: networks/with_blocked_reactions/ of KIMONU_RECORDS in, networks/unblocked/ of KIMONU_OUT out
# (records.py). With another KIMONU_OUT submit with -o <that folder>/networks/unblocked/slurm_%j.log.
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
export PYTHONPATH=$HOME/sd-scip-pricing
PY=$HOME/miniconda3/envs/straindesign/bin/python
export IN=$($PY $BASE/scripts/records.py read WITH_BLOCKED) DEST=$($PY $BASE/scripts/records.py write UNBLOCKED)
ls $IN/*_fixed.crn.json | xargs -P 16 -I{} bash -c '
  f={}; m=$(basename $f _fixed.crn.json)
  BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
  $HOME/miniconda3/envs/straindesign/bin/python $BASE/scripts/write_unblocked_crn.py \
    $f $DEST/${m}_fixed.crn.json --solver scip \
    > $DEST/${m}.log 2>&1 || echo "FAILED $m"
'
echo "ALL DONE $(date -Is)"
