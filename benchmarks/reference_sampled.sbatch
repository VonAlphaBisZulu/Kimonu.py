#!/usr/bin/env bash
#SBATCH --job-name=refsampled
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/reference_sampled/slurm_%j.log
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=$BASE/benchmarks/reference_sampled
U=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm
BUDGET=${BUDGET:-86400}
JOBS=${JOBS:-30}
SEED=${SEED:-1234}
mkdir -p $OUT

# The reference's own workflow: pre_concordant_sampling discards pairs whose activity ratio visibly
# varies across sampled fluxes, then the exact LP test runs on the survivors. Only ever removes
# candidates, so it can lose a true concordant pair -- which is why the exhaustive sweep is run too.
for f in $(ls -S -r $BASE/benchmarks/unblocked/*_fixed.crn.json); do
  basename $f _fixed.crn.json
done | xargs -P $JOBS -I{} bash -c '
  m={}
  BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
  OUT=$BASE/benchmarks/reference_sampled
  export PATH=$HOME/miniconda3/envs/octave/bin:$PATH
  export OCTAVE_HOME=$HOME/miniconda3/envs/octave
  export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
  export UPSTREAM_ALGORITHM=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm
  export REFERENCE_MODE=sampled
  export REFERENCE_SEED='"${SEED:-1234}"'
  export REFERENCE_TIME_BUDGET='"${BUDGET:-86400}"'
  export REFERENCE_MAT_OUT=$OUT/${m}_concordant_sampled.mat
  [ -s "$OUT/reference_${m}.json" ] && { echo "have $m"; exit 0; }
  stdbuf -oL -eL octave --no-gui $BASE/benchmarks/matlab_reference/run_reference.m \
    $BASE/benchmarks/unblocked/${m}_fixed.crn.json $OUT/reference_${m}.json \
    > $OUT/reference_${m}.log 2>&1
  echo "done $m rc=$?"
'
echo "ALL SAMPLED CELLS DONE $(date -Is)"
