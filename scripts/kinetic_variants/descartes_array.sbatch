#!/usr/bin/env bash
#SBATCH --job-name=descartes
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-3
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_logs/descartes_%A_%a.log
# one model per task: Jacobian route vs iterated route, both uptake modes
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models_gap.txt)
/usr/bin/time -v $HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/descartes_vs_tangent.py "$M"
