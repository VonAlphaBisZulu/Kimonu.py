#!/usr/bin/env bash
#SBATCH --job-name=reduce-gap
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts/kinetic_variants/logs/%x_%j.log
# $1 model, $2 target metabolite, $3 output json
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/reduce_gap.py "$@"
