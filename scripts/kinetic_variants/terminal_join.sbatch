#!/usr/bin/env bash
#SBATCH --job-name=termjoin
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-06:00:00
#SBATCH --array=0-32
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_terminal/logs/%x_%A_%a.log
# one model per task, attached uptakes; results as json lines in /scratch/.../kimonu_terminal/
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
mapfile -t MODELS < kinetic_variants/models33.txt
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/terminal_join.py "$M" \
    > /scratch/sadreyev_lab/pschneider/kimonu_terminal/$M.json
