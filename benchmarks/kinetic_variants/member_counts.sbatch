#!/usr/bin/env bash
#SBATCH --job-name=member_counts
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-04:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_member_counts/logs/mc_%j.log
# all 33 models in turn: member counts of the iterated route, and of the tangent route on the gap models
#   sbatch kinetic_variants/member_counts.sbatch [attached]   (MC_OUT, UM_DIR as in member_counts.py)
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
export PYTHONPATH=$PWD:${SD_TREE:-$HOME/sd-kernel-nonzeros} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
bash kinetic_variants/provenance.sh "${MC_OUT:-$PWD/kinetic_variants/member_counts}"
echo "$(hostname) | kimonu $(git -C .. rev-parse --short HEAD) | straindesign $(git -C ${SD_TREE:-$HOME/sd-kernel-nonzeros} rev-parse --short HEAD) | $(date -Is)"
while read -r M; do
  $HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/member_counts.py "$M" "$@" || echo "FAILED $M"
done < kinetic_variants/models33.txt
date -Is
