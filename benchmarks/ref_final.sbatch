#!/usr/bin/env bash
#SBATCH --job-name=reffinal
#SBATCH --nodes=1

#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --partition=sadreyevlab
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/reference_final/full_%j.log
set -u
# the log directory must exist before Slurm opens the log; see FINDINGS
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=$BASE/benchmarks/reference_final
mkdir -p $OUT
export PATH=$HOME/miniconda3/envs/octave/bin:$PATH
export OCTAVE_HOME=$HOME/miniconda3/envs/octave
export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
export UPSTREAM_ALGORITHM=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm  # fixed concordance (4 commits) -- NOT the pinned worktree
export REFERENCE_MODE=exhaustive
# no REFERENCE_TIME_BUDGET: run to completion, to widen the correctness ground truth beyond
# the four small iIS312 variants
export REFERENCE_MAT_OUT=$OUT/${MODEL}_concordant_fixed.mat
[ -s "$OUT/reference_${MODEL}.json" ] && { echo "have $MODEL"; exit 0; }
stdbuf -oL -eL octave --no-gui $BASE/benchmarks/matlab_reference/run_reference.m \
  $BASE/benchmarks/unblocked/${MODEL}_fixed.crn.json $OUT/reference_${MODEL}.json
echo "done $MODEL rc=$?"
