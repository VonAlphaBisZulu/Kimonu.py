#!/usr/bin/env bash
#SBATCH --job-name=tp-local
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-06:00:00
#SBATCH --array=1-51
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_twopoint_local/logs/%x_%A_%a.log
# one relation per task: the sign certificate on a growing local subsystem (two_point_local.py)
# task list: /scratch/sadreyev_lab/pschneider/kimonu_twopoint_local/tasks.txt, lines "<model> <relation index>"
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
read M I < <(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/sadreyev_lab/pschneider/kimonu_twopoint_local/tasks.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/two_point_local.py "$M" "$I" 900
