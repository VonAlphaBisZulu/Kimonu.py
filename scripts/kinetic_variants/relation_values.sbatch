#!/usr/bin/env bash
#SBATCH --job-name=relation_values
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-12:00:00
#SBATCH --array=1-33
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_values/logs/rv_%A_%a.log
# one model per task: exact relation values at given k (round trip, prime k), attached and detached
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/relation_values.py "$M"
