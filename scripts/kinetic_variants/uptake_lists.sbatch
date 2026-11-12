#!/usr/bin/env bash
#SBATCH --job-name=uptake_lists
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --array=1-33
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_uptake_modes/logs/ul_%A_%a.log
# one model per task: iterated ACR / ACRR lists, attached and detached uptakes. Writes results/kimonu/acr_lists_by_uptake_mode/ of KIMONU_OUT (records.py).
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
SD=${SD_TREE:-$HOME/sd-kernel-nonzeros}
export PYTHONPATH=$PWD:$SD OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
OUT=$($HOME/miniconda3/envs/straindesign/bin/python records.py write UPTAKE_LISTS)
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
echo "$(hostname) | $M | kimonu $(git -C .. rev-parse --short HEAD) | straindesign $(git -C $SD describe --tags --always) | $(date -Is)"
bash kinetic_variants/provenance.sh "$OUT"
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/uptake_lists.py "$M" "$OUT"
