#!/usr/bin/env bash
#SBATCH --job-name=kimonu_preprocess
#SBATCH --partition=defq
#SBATCH --exclude=node[06-07,09-16]
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_timing_seeds/fva/preprocess_%j.log
# What a user pays before the analysis on the two largest networks: load the SBML, split reversible reactions, find
# the blocked reactions with one LP, build N = Y A (about 12 s on iAF1260b and iS_1188, one core, exclusive node).
# Ran 2026-11-06 (job 570671) with StrainDesign 1.19.2.
#   sbatch kinetic_variants/preprocess_time.sbatch        # from benchmarks/
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-kernel-nonzeros:../src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=$HOME/miniconda3/envs/straindesign/bin/python; O=${FVA_OUT:-/scratch/sadreyev_lab/pschneider/kimonu_timing_seeds/fva}
mkdir -p $O
CORE=$($PY -c "import os; print(min(os.sched_getaffinity(0)))"); hostname
for m in iAF1260b iS_1188; do
  [ -s $O/$m.xml ] || $PY kinetic_variants/crn_to_sbml.py sweep/${m}_fixed.crn.json $O/$m.xml 2>&1 | tail -1
  for r in 1 2 3; do /usr/bin/time -f "peak %M kB, wall %e s" taskset -c $CORE $PY kinetic_variants/preprocess_time.py $O/$m.xml 2>&1 | grep -v -i "warn\|^$" | tail -2 | cut -c1-260; done
done
