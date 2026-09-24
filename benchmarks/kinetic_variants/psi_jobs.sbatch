#!/usr/bin/env bash
#SBATCH --job-name=psi-kin
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/kinetic_variants/logs/%x_%j.log
# $1 = script (psi_nullspace.py | psi_validate.py), rest = its arguments
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
s=$1; shift
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/$s "$@"
