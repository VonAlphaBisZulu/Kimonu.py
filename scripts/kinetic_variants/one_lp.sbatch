#!/usr/bin/env bash
#SBATCH --job-name=kimonu_one_lp
#SBATCH --partition=defq
#SBATCH --exclude=node[06-07,09-16]
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_timing_seeds/fva/one_lp_%j.log
# One LP for blocked reactions, against the FVA: one core, one thread, nothing else on the node. The result lines go to
# results/kimonu/one_lp_vs_fva/one_lp_results.log of KIMONU_OUT (records.py), the networks come from KIMONU_RECORDS.
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=$HOME/miniconda3/envs/straindesign/bin/python
O=$($PY records.py write ONE_LP); W=$($PY records.py read WITH_BLOCKED); U=$($PY records.py read UNBLOCKED)
CORE=$($PY -c "import os; print(min(os.sched_getaffinity(0)))")
for f in $W/*_fixed.crn.json; do
  m=$(basename $f _fixed.crn.json)
  taskset -c $CORE $PY kinetic_variants/one_lp.py $f $U/${m}_fixed.crn.json 2>&1 | tail -1 | cut -c1-400
done > $O/one_lp_results.log
