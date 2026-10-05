#!/usr/bin/env bash
#SBATCH --job-name=cocoa-loop
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_cocoa_loop/logs/%x_%A_%a.log
# sbatch --array=0-32 cocoa_loop.sbatch : Kimonu concordance + iterated route, then COCOA's kinetic stage on it
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=/scratch/sadreyev_lab/pschneider/kimonu_cocoa_loop
mapfile -t MODELS < $BASE/benchmarks/cocoa_loop/models.txt
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
cd $BASE
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 JULIA_NUM_THREADS=8
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia PATH=$HOME/.local/bin:$PATH
[ -s $OUT/kimonu_$M.json ] || $HOME/miniconda3/envs/straindesign/bin/python -u benchmarks/cocoa_loop/kimonu_side.py $M $OUT/kimonu_$M.json
[ -s $OUT/cocoa_$M.json ] || julia --project=benchmarks/julia_v110 benchmarks/cocoa_loop/cocoa_side.jl \
    benchmarks/unblocked/${M}_fixed.crn.json $OUT/kimonu_$M.json $OUT/cocoa_$M.json
echo "--- $M done $(date -Is)"
