#!/usr/bin/env bash
#SBATCH --job-name=reffinal-s
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# peak RAM is the sampling prefilter's m x 10 x m ratio array (~7 * complexes^2 * 80 B);
# override --mem per model on the sbatch command line
#SBATCH --partition=sadreyevlab
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/reference_final_sampled/s_%A_%a.log
set -u
# the log directory must exist before Slurm opens the log; see FINDINGS
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=$BASE/benchmarks/reference_final_sampled; mkdir -p $OUT
mapfile -t CELLS < $BASE/benchmarks/final_sampled_cells.txt
spec=${CELLS[$SLURM_ARRAY_TASK_ID]}; m=${spec%%:*}; S=${spec##*:}
export PATH=$HOME/miniconda3/envs/octave/bin:$PATH OCTAVE_HOME=$HOME/miniconda3/envs/octave
export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
export UPSTREAM_ALGORITHM=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm   # fixed concordance -- NOT the pinned worktree
export REFERENCE_MODE=sampled REFERENCE_SEED=$S
export REFERENCE_MAT_OUT=$OUT/${m}_seed${S}.mat
[ -s "$OUT/reference_${m}_seed${S}.json" ] && { echo "have $m $S"; exit 0; }
stdbuf -oL -eL octave --no-gui $BASE/benchmarks/matlab_reference/run_reference.m \
  $BASE/benchmarks/unblocked/${m}_fixed.crn.json $OUT/reference_${m}_seed${S}.json
echo "done $m $S rc=$?"
