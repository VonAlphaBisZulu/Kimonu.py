#!/usr/bin/env bash
#SBATCH --job-name=refseed
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/reference_sampled_seeds/rs_%A_%a.log
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=$BASE/benchmarks/reference_sampled_seeds
mkdir -p $OUT
mapfile -t CELLS < $BASE/benchmarks/refseed_cells.txt
spec=${CELLS[$SLURM_ARRAY_TASK_ID]}
m=${spec%%:*}; S=${spec##*:}
export PATH=$HOME/miniconda3/envs/octave/bin:$PATH
export OCTAVE_HOME=$HOME/miniconda3/envs/octave
export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
export UPSTREAM_ALGORITHM=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm
export REFERENCE_MODE=sampled
export REFERENCE_SEED=$S
export REFERENCE_TIME_BUDGET=86400
export REFERENCE_MAT_OUT=$OUT/${m}_seed${S}.mat
[ -s "$OUT/reference_${m}_seed${S}.json" ] && { echo "have $m $S"; exit 0; }
stdbuf -oL -eL octave --no-gui $BASE/benchmarks/matlab_reference/run_reference.m \
  $BASE/benchmarks/unblocked/${m}_fixed.crn.json $OUT/reference_${m}_seed${S}.json
echo "done $m seed $S rc=$?"
