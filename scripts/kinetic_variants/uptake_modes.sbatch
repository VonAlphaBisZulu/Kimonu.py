#!/usr/bin/env bash
#SBATCH --job-name=uptake_modes
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-08:00:00
#SBATCH --array=1-33
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_uptake_modes/logs/um_%A_%a.log
# one model per task: both routes, gap relations and g_rel, attached and detached uptakes. UM_OUT replaces the output
# directory kinetic_variants/uptake_modes, SD_TREE the StrainDesign tree (default v1.19.2).
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
SD=${SD_TREE:-$HOME/sd-kernel-nonzeros}
export PYTHONPATH=$PWD:$SD OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
echo "$(hostname) | $M | kimonu $(git -C .. rev-parse --short HEAD) | straindesign $(git -C $SD describe --tags --always) | $(date -Is)"
bash kinetic_variants/provenance.sh "${UM_OUT:-$PWD/kinetic_variants/uptake_modes}"
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/uptake_modes.py "$M"
