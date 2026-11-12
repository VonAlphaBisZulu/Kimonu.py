#!/usr/bin/env bash
#SBATCH --job-name=injectivity
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-05:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_injectivity/logs/%x_%j.log
# sbatch injectivity.sbatch <model> dims|inject|quotient [time_limit_s] [threads]
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/injectivity_gap.py "$@"
