#!/usr/bin/env bash
#SBATCH --job-name=lacr-ladder
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=0-03:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_local_acr/logs/%x_%j.log
# sbatch local_acr_ladder.sbatch <network> [limit_s]: symbolic local-ACR criterion vs exact nullspace, one network
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/local_acr_symbolic.py ladder "$1" "${2:-3600}"
