#!/usr/bin/env bash
#SBATCH --job-name=gap_roots
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-04:00:00
#SBATCH --array=1-15
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_gaproots/logs/gr_%A_%a.log
# The open relations on the three gap networks that follow from a few roots (Results, "at most 5, 1 and 2 open
# questions"). One case per task from gap_roots_cases.txt: DETACH flag, model, roots (metabolite = ACR, a/b = ACRR,
# "none" = no root). Ran 2026-11-02 (job 568497) with StrainDesign ~/sd-scip-pricing; the
# default tree is now StrainDesign 1.19.2 ($HOME/sd-kernel-nonzeros).
#   sbatch kinetic_variants/gap_roots.sbatch        # from scripts/
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
SD=${SD_TREE:-$HOME/sd-kernel-nonzeros}
export PYTHONPATH=$SD OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
OUT=${GAP_OUT:-/scratch/sadreyev_lab/pschneider/kimonu_gaproots}
mkdir -p $OUT
echo "$(hostname) | kimonu $(git -C .. describe --always --dirty) | straindesign $(git -C $SD rev-parse --short HEAD) | $(date -Is)" >&2
bash kinetic_variants/provenance.sh "$OUT"
read -r D M ROOTS <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/gap_roots_cases.txt)"
DETACH=$D $HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/gap_roots.py "$M" $ROOTS \
  > $OUT/case_${SLURM_ARRAY_TASK_ID}.json
