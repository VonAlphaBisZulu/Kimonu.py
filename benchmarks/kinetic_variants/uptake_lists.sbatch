#!/usr/bin/env bash
#SBATCH --job-name=uptake_lists
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --array=1-33
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_uptake_modes/logs/ul_%A_%a.log
# one model per task: iterated ACR / ACRR lists, attached and detached uptakes
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/uptake_lists.py "$M" /scratch/sadreyev_lab/pschneider/kimonu_uptake_modes/lists
