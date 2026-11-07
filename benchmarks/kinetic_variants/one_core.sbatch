#!/usr/bin/env bash
#SBATCH --job-name=kimonu_1cpu
#SBATCH --partition=defq
#SBATCH --exclude=node[06-09,16]
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_lean/logs/one_%A_%a.log
# one model and replicate per task on one reserved core and 8 GB: does Kimonu need more?
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
L=/scratch/sadreyev_lab/pschneider/kimonu_lean
read A M R < <(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/one_core_cells.txt)
PY=$HOME/miniconda3/envs/straindesign/bin/python
echo "$A $M rep$R node=$(hostname) cpus=$(nproc)"
if [ "$A" = det ]; then
  /usr/bin/time -v $PY run_kimonu_crn.py unblocked/${M}_fixed.crn.json $L/det/kimonu_${M}_rep${R}.json --repeats 1 --assume-unblocked --solver scip --detach-uptakes > $L/det/kimonu_${M}_rep${R}.log 2>&1
else
  /usr/bin/time -v $PY run_kimonu_refmode.py $M $L/ref/kimonu_${M}_rep${R}.json > $L/ref/kimonu_${M}_rep${R}.log 2>&1
fi
echo "rc=$?"
