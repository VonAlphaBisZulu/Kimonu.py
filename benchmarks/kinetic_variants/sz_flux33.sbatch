#!/usr/bin/env bash
#SBATCH --job-name=sz_flux33
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-33
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_sz_flux_33/logs/sz_%A_%a.log
# All 33 networks: both routes with the default r* and with the uniformly drawn one (2^64), StrainDesign v1.19.2.
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-kernel-nonzeros OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SZ_OUT=/scratch/sadreyev_lab/pschneider/kimonu_sz_flux_33
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/sz_flux.py "$M" 64
