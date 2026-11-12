#!/usr/bin/env bash
#SBATCH --job-name=seed_inv
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-04:00:00
#SBATCH --array=1-33
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_seeds/logs/s_%A_%a.log
# Kimonu's results under three seeds of the reference flux, one network per task (fig. S4, "three seeds identical").
# Writes $SEED_OUT/<model>.json (default /scratch/sadreyev_lab/pschneider/kimonu_seeds). Ran 2026-11-02 (job 568581)
# with StrainDesign ~/sd-scip-pricing; set SD_TREE=$HOME/sd-kernel-nonzeros for StrainDesign 1.19.2.
#   sbatch kinetic_variants/seed_invariance.sbatch        # from scripts/
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
SD=${SD_TREE:-$HOME/sd-scip-pricing}
export PYTHONPATH=$SD OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
echo "$(hostname) | kimonu $(git -C .. describe --always --dirty) | straindesign $(git -C $SD rev-parse --short HEAD) | $(date -Is)"
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models33.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/seed_invariance.py "$M"
