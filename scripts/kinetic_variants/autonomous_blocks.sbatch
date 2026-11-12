#!/usr/bin/env bash
#SBATCH --job-name=autoblocks
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-12:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_autoblocks/logs/%x_%A_%a.log
# sbatch --array=0-32 autonomous_blocks.sbatch [extra args]; one model per task from models33.txt
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
mapfile -t MODELS < kinetic_variants/models33.txt
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/autonomous_blocks.py "$M" "$@" \
    > /scratch/sadreyev_lab/pschneider/kimonu_autoblocks/$M.json
