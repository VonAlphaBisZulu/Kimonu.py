#!/usr/bin/env bash
#SBATCH --job-name=second_order
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-06:00:00
#SBATCH --array=1-3
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_second_order/logs/so_%A_%a.log
# one gap model per task: second-order test of the tangent-only relations
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models_gap.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/${SO_SCRIPT:-second_order.py} "$M"
