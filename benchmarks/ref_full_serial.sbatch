#!/usr/bin/env bash
#SBATCH --job-name=reffull
#SBATCH --nodes=1

#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=sadreyevlab
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/reference_full/full_%j.log
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=$BASE/benchmarks/reference_full
mkdir -p $OUT
export PATH=$HOME/miniconda3/envs/octave/bin:$PATH
export OCTAVE_HOME=$HOME/miniconda3/envs/octave
export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
# Pinned to the worktree at 088d827, the code that produced every completed reference run.
# The main checkout carries the proposed upstream fixes and must not be used for benchmarks.
export UPSTREAM_ALGORITHM=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm_benchmark
export REFERENCE_MODE=exhaustive
# no REFERENCE_TIME_BUDGET: run to completion, to widen the correctness ground truth beyond
# the four small iIS312 variants
export REFERENCE_MAT_OUT=$OUT/${MODEL}_concordant_fixed.mat
[ -s "$OUT/reference_${MODEL}.json" ] && { echo "have $MODEL"; exit 0; }
stdbuf -oL -eL octave --no-gui $BASE/benchmarks/matlab_reference/run_reference.m \
  $BASE/benchmarks/unblocked/${MODEL}_fixed.crn.json $OUT/reference_${MODEL}.json
echo "done $MODEL rc=$?"
