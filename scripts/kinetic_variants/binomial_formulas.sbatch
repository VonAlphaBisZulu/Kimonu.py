#!/usr/bin/env bash
#SBATCH --job-name=binomial_formulas
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-12:00:00
#SBATCH --array=1-33
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_binomial/logs/bf_%A_%a.log
# one model per task: monomial formulas from flux coupling, checked against exact values
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/binomial_formulas.py "$M"
