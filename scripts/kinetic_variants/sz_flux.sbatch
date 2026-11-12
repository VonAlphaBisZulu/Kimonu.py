#!/usr/bin/env bash
#SBATCH --job-name=sz_flux
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-8
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_sz_flux/logs/sz_%A_%a.log
# one model per task (sz_sample.txt): both routes with the default r* and with a uniformly drawn one (2^64)
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=${SD_TREE:-$HOME/sd-kernel-nonzeros} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/sz_sample.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/sz_flux.py "$M" 64
