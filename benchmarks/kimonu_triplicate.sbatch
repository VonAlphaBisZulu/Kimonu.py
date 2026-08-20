#!/usr/bin/env bash
#SBATCH --job-name=kimtrip
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/triplicate/slurm_%j.log
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=$BASE/benchmarks/${OUTDIR:-triplicate}
CRN=$BASE/benchmarks/unblocked
PY=$HOME/miniconda3/envs/straindesign/bin/python
export PYTHONPATH=${SD_BRANCH:-$HOME/sd-scip-pricing}
SOLVER=${SOLVER:-scip}
REPS=${REPS:-3}
JOBS=${JOBS:-24}
mkdir -p $OUT

# One cell per (model, replicate), each single-core, run through xargs so the node is filled
# without any cell seeing more than one CPU. Replicates are separate processes rather than
# --repeats inside one, so each carries its own interpreter and solver start-up like a real run.
for crn in $CRN/*_fixed.crn.json; do
  m=$(basename $crn _fixed.crn.json)
  for r in $(seq 1 $REPS); do
    echo "$m $r"
  done
done | xargs -P $JOBS -n 2 bash -c '
  m=$0; r=$1
  BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
  OUT=$BASE/benchmarks/${OUTDIR:-triplicate}
  /usr/bin/time -v $HOME/miniconda3/envs/straindesign/bin/python \
    $BASE/benchmarks/run_kimonu_crn.py $BASE/benchmarks/unblocked/${m}_fixed.crn.json \
    $OUT/kimonu_${m}_rep${r}.json --repeats 1 --assume-unblocked --solver '"${SOLVER:-scip}"' \
    > $OUT/kimonu_${m}_rep${r}.log 2>&1
  echo "done $m rep$r rc=$?"
'
echo "ALL CELLS DONE $(date -Is)"
