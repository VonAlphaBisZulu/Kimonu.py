#!/usr/bin/env bash
#SBATCH --job-name=cert_modp
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-33
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_cert_modp/logs/cm_%A_%a.log
# one model per task: both routes again at a uniform random point mod 2^61 - 1, against the exact results
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=${SD_TREE:-$HOME/sd-kernel-nonzeros} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/cert_modp.py "$M" 0
