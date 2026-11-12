#!/usr/bin/env bash
#SBATCH --job-name=kimonu_fva_time
#SBATCH --partition=defq
#SBATCH --exclude=node[06-07,09-16]
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_timing_seeds/fva/fva_%j.log
# Wall time of the blocked-reaction FVA that write_unblocked_crn.py runs (one core, nothing else on the node), to set
# against the one LP of one_lp.sbatch. Ran 2026-11-06 (job 570669) with StrainDesign 1.19.2.
#   sbatch kinetic_variants/fva_time.sbatch        # from scripts/
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-kernel-nonzeros OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=$HOME/miniconda3/envs/straindesign/bin/python
O=$($PY records.py write EXPERIMENTS fva_time); W=$($PY records.py read WITH_BLOCKED)   # records.py
CORE=$($PY -c "import os; print(min(os.sched_getaffinity(0)))")
echo "$(hostname) core $CORE"
for f in $W/*_fixed.crn.json; do
  m=$(basename $f _fixed.crn.json)
  /usr/bin/time -f "$m %e s %M kB" taskset -c $CORE $PY kinetic_variants/fva_one_process.py write_unblocked_crn.py $f $O/$m.crn.json --solver scip 2>&1 | tail -2 | cut -c1-200
done
