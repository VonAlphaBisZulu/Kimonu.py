#!/usr/bin/env bash
#SBATCH --job-name=ref-finish
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_records/results/reference/exhaustive/finish_%x_%j.log
# Finish an exhaustive reference run: re-test the pairs it rejected (scaling-bound fix) when RETEST=1,
# stage the .mat for the R pipeline, run the R stages (kinetic modules, ACR/ACRR) and dump the modules.
# sbatch --export=ALL,MODEL=<model>,RETEST=0|1 ref_finish.sbatch
# Reads results/reference/exhaustive/ of KIMONU_RECORDS and writes there in KIMONU_OUT (records.py); the log line names
# the records folder of the cluster, so with another KIMONU_OUT submit with -o <that folder>/finish_%x_%j.log.
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
PY=$HOME/miniconda3/envs/straindesign/bin/python
export PATH=$HOME/miniconda3/envs/octave/bin:$PATH OCTAVE_HOME=$HOME/miniconda3/envs/octave
export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
export UPSTREAM_ALGORITHM=${UPSTREAM_ALGORITHM:-$BASE/reference/Upstream_Algorithm}   # bash reference/fetch_upstream.sh
export PYTHONPATH=$HOME/sd-scip-pricing
IN=$($PY $BASE/records.py read REF_EXHAUSTIVE); OUT=$($PY $BASE/records.py write REF_EXHAUSTIVE)
mkdir -p $OUT/retest $OUT/acr/Results/concordant
stem=${MODEL}_concordant_fixed
if [ "${RETEST:-0}" = 1 ]; then
  $PY $BASE/matlab_reference/retest_pairs.py $MODEL
  (cd $BASE/matlab_reference && stdbuf -oL octave --no-gui retest_rejected.m $MODEL) > $OUT/retest/retest_$MODEL.log 2>&1
  cp $OUT/retest/${stem}_retest.mat $OUT/acr/Results/concordant/
  stem=${stem}_retest
else
  cp $IN/${stem}.mat $OUT/acr/Results/concordant/
fi
bash $BASE/reference/run_chain.sh exhaustive $stem
echo "finished $MODEL ($stem) $(date -Is)"
