#!/usr/bin/env bash
#SBATCH --job-name=float_exactcmp
#SBATCH --partition=defq
#SBATCH --exclude=node[01-09,16]
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_float_four/logs/kx_%A_%a.log
# one model and part per task (list file as $1): exact kernel, exact tolerance comparison; detached uptakes
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=${SD_TREE:-$HOME/sd-kernel-nonzeros} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
read M P < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/float_exactcmp.py "$M" "$P"
