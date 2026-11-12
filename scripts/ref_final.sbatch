#!/usr/bin/env bash
#SBATCH --job-name=reffinal
#SBATCH --nodes=1

#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --partition=sadreyevlab
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_records/results/reference/exhaustive/full_%j.log
# Frozen harness of the paper (the reference of Langary et al. 2025). Outputs go to results/reference/exhaustive/ of KIMONU_OUT (records.py);
# the log line names the records folder of the cluster, so with another KIMONU_OUT submit with -o <that folder>/full_%j.log.
set -u
# the log directory must exist before Slurm opens the log; see FINDINGS
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
PY=$HOME/miniconda3/envs/straindesign/bin/python
OUT=$($PY $BASE/scripts/records.py write REF_EXHAUSTIVE)
NET=$($PY $BASE/scripts/records.py read UNBLOCKED)
export PATH=$HOME/miniconda3/envs/octave/bin:$PATH
export OCTAVE_HOME=$HOME/miniconda3/envs/octave
export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
export UPSTREAM_ALGORITHM=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm  # fixed concordance (4 commits) -- NOT the pinned worktree
export REFERENCE_MODE=exhaustive
# no REFERENCE_TIME_BUDGET: run to completion, to widen the correctness ground truth beyond
# the four small iIS312 variants
export REFERENCE_MAT_OUT=$OUT/${MODEL}_concordant_fixed.mat
[ -s "$OUT/reference_${MODEL}.json" ] && { echo "have $MODEL"; exit 0; }
stdbuf -oL -eL octave --no-gui $BASE/scripts/matlab_reference/run_reference.m \
  $NET/${MODEL}_fixed.crn.json $OUT/reference_${MODEL}.json
echo "done $MODEL rc=$?"
