#!/usr/bin/env bash
#SBATCH --job-name=second_order_modp
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-06:00:00
#SBATCH --array=1-3
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_second_order/logs/so_%A_%a.log
# Exact second-order test of the relations left open on the three gap networks, over GF(2^31 - 1) (Methods,
# Supplementary Text). One gap network per task; writes $SO_OUT/<model>_modp.json. Ran 2026-10-28 (job 567790)
# with StrainDesign ~/sd-scip-pricing; the default tree is now StrainDesign 1.19.2
# ($HOME/sd-kernel-nonzeros).
#   sbatch kinetic_variants/second_order_modp.sbatch        # from benchmarks/
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
SD=${SD_TREE:-$HOME/sd-kernel-nonzeros}
export PYTHONPATH=$SD OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
echo "$(hostname) | kimonu $(git -C .. describe --always --dirty) | straindesign $(git -C $SD rev-parse --short HEAD) | $(date -Is)" >&2
bash kinetic_variants/provenance.sh "${SO_OUT:-/scratch/sadreyev_lab/pschneider/kimonu_second_order}"
M=$(sed -n "${SLURM_ARRAY_TASK_ID}p" kinetic_variants/models_gap.txt)
$HOME/miniconda3/envs/straindesign/bin/python -u kinetic_variants/second_order_modp.py "$M"
