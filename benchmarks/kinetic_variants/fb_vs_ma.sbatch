#!/usr/bin/env bash
#SBATCH --job-name=fb_vs_ma
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-06:00:00
#SBATCH --array=1-33%16
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_fb_vs_ma/logs/fm_%A_%a.log
# one model per task: flux coupling and concordance against their mass-action counterparts, both uptake modes
#   sbatch kinetic_variants/fb_vs_ma.sbatch   (FB_OUT, UM_DIR as in fb_vs_ma.py)
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$PWD:${SD_TREE:-$HOME/sd-kernel-nonzeros} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
bash kinetic_variants/provenance.sh "${FB_OUT:-$PWD/kinetic_variants/fb_vs_ma}"
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
echo "$(hostname) | $M | kimonu $(git -C .. rev-parse --short HEAD) | straindesign $(git -C ${SD_TREE:-$HOME/sd-kernel-nonzeros} rev-parse --short HEAD) | $(date -Is)"
for mode in detached attached; do
  $HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/fb_vs_ma.py "$M" $mode || echo "FAILED $M $mode"
done
date -Is
