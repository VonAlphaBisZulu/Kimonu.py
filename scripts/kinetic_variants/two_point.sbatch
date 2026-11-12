#!/usr/bin/env bash
#SBATCH --job-name=two-point
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-06:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_twopoint/logs/%x_%A_%a.log
# prepare:  sbatch two_point.sbatch <model> prepare
# relation: sbatch --array=0-K two_point.sbatch <model> relation [search s] [milp s]   (index = array task)
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
M=$1; CMD=$2; shift 2
if [ "$CMD" = relation ]; then set -- "$SLURM_ARRAY_TASK_ID" "$@"; fi
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/two_point_gap.py "$M" "$CMD" "$@"
