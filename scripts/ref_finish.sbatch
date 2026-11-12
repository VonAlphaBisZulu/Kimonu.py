#!/usr/bin/env bash
#SBATCH --job-name=ref-finish
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts/reference_final/finish_%x_%j.log
# Finish an exhaustive reference run: re-test the pairs it rejected (scaling-bound fix) when RETEST=1,
# stage the .mat for the R pipeline, run the R stages (kinetic modules, ACR/ACRR) and dump the modules.
# sbatch --export=ALL,MODEL=<model>,RETEST=0|1 ref_finish.sbatch
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PATH=$HOME/miniconda3/envs/octave/bin:$PATH OCTAVE_HOME=$HOME/miniconda3/envs/octave
export LP_SERVER_PYTHON=$HOME/miniconda3/envs/straindesign/bin/python3
export UPSTREAM_ALGORITHM=$HOME/pipelines/2026-08_Kimonu/Upstream_Algorithm
export PYTHONPATH=$HOME/sd-scip-pricing
stem=${MODEL}_concordant_fixed
if [ "${RETEST:-0}" = 1 ]; then
  $HOME/miniconda3/envs/straindesign/bin/python $BASE/matlab_reference/retest_pairs.py $MODEL
  (cd $BASE/matlab_reference && stdbuf -oL octave --no-gui retest_rejected.m $MODEL) > $BASE/reference_final/retest/retest_$MODEL.log 2>&1
  cp $BASE/reference_final/retest/${stem}_retest.mat $BASE/reference_final/acr/Results/concordant/
  stem=${stem}_retest
else
  cp $BASE/reference_final/${stem}.mat $BASE/reference_final/acr/Results/concordant/
fi
bash $BASE/reference_final/acr/run_chain.sh $stem
echo "finished $MODEL ($stem) $(date -Is)"
