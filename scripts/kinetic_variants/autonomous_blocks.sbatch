#!/usr/bin/env bash
#SBATCH --job-name=autoblocks
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-12:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_autoblocks/logs/%x_%A_%a.log
# sbatch --array=0-32 autonomous_blocks.sbatch [extra args]; one model per task from models33.txt, one JSON line per model in
# results/kimonu/gap_experiments/autonomous_blocks/ of KIMONU_OUT (records.py); autonomous_blocks_33.jsonl there is the 33 lines
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
export PYTHONPATH=$HOME/sd-scip-pricing OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
mapfile -t MODELS < kinetic_variants/models33.txt
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/autonomous_blocks.py "$M" "$@" \
    > "$($HOME/miniconda3/envs/straindesign/bin/python records.py write AUTONOMOUS_BLOCKS)/$M.json"
