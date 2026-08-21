#!/usr/bin/env bash
#SBATCH --job-name=refsweep
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/reference/slurm_%j.log
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=$BASE/benchmarks/reference
U=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm
BUDGET=${BUDGET:-86400}
JOBS=${JOBS:-30}
mkdir -p $OUT

# The reference is serial by construction, so the node is filled by running many models at once
# rather than by giving any one model more cores. Ordered smallest first so results accumulate.
# Each model's linprog calls go to its own persistent CPLEX worker; glpk both stalled indefinitely
# on some complexes and returned spurious zeros, inflating the balanced set by 26% on iIT341.
for f in $(ls -S -r $BASE/benchmarks/unblocked/*_fixed.crn.json); do
  basename $f _fixed.crn.json
done | xargs -P $JOBS -I{} bash -c '
  m={}
  BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
  OUT=$BASE/benchmarks/reference
  U=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm
  export PATH=$HOME/miniconda3/envs/octave/bin:$PATH
  export OCTAVE_HOME=$HOME/miniconda3/envs/octave
  export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
  export UPSTREAM_ALGORITHM=$U
  export REFERENCE_MODE=exhaustive
  export REFERENCE_TIME_BUDGET='"${BUDGET:-86400}"'
  export REFERENCE_MAT_OUT=$OUT/${m}_concordant_fixed.mat
  stdbuf -oL -eL octave --no-gui $BASE/benchmarks/matlab_reference/run_reference.m \
    $BASE/benchmarks/unblocked/${m}_fixed.crn.json $OUT/reference_${m}.json \
    > $OUT/reference_${m}.log 2>&1
  echo "done $m rc=$?"
'
echo "ALL REFERENCE CELLS DONE $(date -Is)"
