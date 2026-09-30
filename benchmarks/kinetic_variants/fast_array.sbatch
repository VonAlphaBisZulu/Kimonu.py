#!/usr/bin/env bash
#SBATCH --job-name=tangent-fast
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-33
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/kinetic_variants/logs/fast_%A_%a.log
# one model per task: Jacobian route vs iterated route, both uptake modes
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
/usr/bin/time -v $HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/fast_vs_old.py "$M"
