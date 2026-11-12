#!/usr/bin/env bash
#SBATCH --job-name=float_four
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-297
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_float_four/logs/kf_%A_%a.log
# one model and part per task (float_four_parts.txt): the tangent arms, or a quarter of the tolerance grid of one iterated arm; detached uptakes
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
read M P < <(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/float_four_parts.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/float_four.py "$M" "$P"
