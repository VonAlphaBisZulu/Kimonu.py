#!/usr/bin/env bash
#SBATCH --job-name=reffinal-s
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# peak RAM is the sampling prefilter's m x 10 x m ratio array (~7 * complexes^2 * 80 B);
# override --mem per model on the sbatch command line
#SBATCH --partition=sadreyevlab
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_records/results/reference/sampled/s_%A_%a.log
# Frozen harness of the paper (the reference of Langary et al. 2025). Outputs go to results/reference/sampled/ of KIMONU_OUT (records.py);
# the log line names the records folder of the cluster, so with another KIMONU_OUT submit with -o <that folder>/s_%A_%a.log.
set -u
# the log directory must exist before Slurm opens the log; see FINDINGS
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
PY=$HOME/miniconda3/envs/straindesign/bin/python
OUT=$($PY $BASE/scripts/records.py write REF_SAMPLED)
NET=$($PY $BASE/scripts/records.py read UNBLOCKED)
mapfile -t CELLS < $BASE/scripts/final_sampled_cells.txt
spec=${CELLS[$SLURM_ARRAY_TASK_ID]}; m=${spec%%:*}; S=${spec##*:}
export PATH=$HOME/miniconda3/envs/octave/bin:$PATH OCTAVE_HOME=$HOME/miniconda3/envs/octave
export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
export UPSTREAM_ALGORITHM=${UPSTREAM_ALGORITHM:-$BASE/scripts/reference/Upstream_Algorithm}   # the fork at 8a6ed17, bash scripts/reference/fetch_upstream.sh
export REFERENCE_MODE=sampled REFERENCE_SEED=$S
export REFERENCE_MAT_OUT=$OUT/${m}_seed${S}.mat
[ -s "$OUT/reference_${m}_seed${S}.json" ] && { echo "have $m $S"; exit 0; }
stdbuf -oL -eL octave --no-gui $BASE/scripts/matlab_reference/run_reference.m \
  $NET/${m}_fixed.crn.json $OUT/reference_${m}_seed${S}.json
echo "done $m $S rc=$?"
