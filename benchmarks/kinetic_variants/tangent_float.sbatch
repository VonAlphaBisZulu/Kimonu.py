#!/usr/bin/env bash
#SBATCH --job-name=tangent-float
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-06:00:00
#SBATCH --array=1-33%11
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_tangent_float/logs/%x_%A_%a.log
# exact vs float64 tangent route, one model per task
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/tangent_float.py "$M"
