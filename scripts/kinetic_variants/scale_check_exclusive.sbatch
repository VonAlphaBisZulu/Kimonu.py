#!/usr/bin/env bash
#SBATCH --job-name=kr_scale_excl
#SBATCH --partition=sadreyevlab
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/logs/%x_%j.log
# S8: LP point rounded at 10^6 against scale 100 and unrounded, iNJ661, three runs, one process pinned to one core.
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$PWD:$HOME/sd-kernel-nonzeros OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=$HOME/miniconda3/envs/straindesign/bin/python
OUT=$($PY records.py write LP_ROUNDING exclusive)   # records.py
CORE=$($PY -c "import os; print(min(os.sched_getaffinity(0)))")
echo "$(hostname) core $CORE | kimonu $(git -C .. rev-parse --short HEAD) | $(date -Is)"
bash kinetic_variants/provenance.sh $OUT hw
for R in 1 2 3; do
  /usr/bin/time -v taskset -c $CORE $PY kinetic_variants/scale_check.py iNJ661 $OUT/iNJ661_rep$R.json > $OUT/iNJ661_rep$R.log 2>&1
  echo "rep $R rc=$?"; grep -E 'tangent|iterated' $OUT/iNJ661_rep$R.log | grep ' s '
done
