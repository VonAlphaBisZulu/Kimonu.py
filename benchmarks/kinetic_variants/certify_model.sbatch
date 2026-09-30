#!/usr/bin/env bash
#SBATCH --job-name=certify
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/kinetic_variants/logs/%x_%j.log
# $1 model, $2 metabolite, $3 time limit
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/certify_model.py "$@"
