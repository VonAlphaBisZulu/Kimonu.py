#!/usr/bin/env bash
#SBATCH --job-name=refmode-mod
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-02:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_trip_2026-11-02/logs/refmod_%j.log
# Kinetic modules of the reference-equivalent mode on all 33 networks, six networks at a time, for the coupled pairs
# of Figure 4C, D (figures/pair_sets.py reads $OUT through KIMONU_REFMOD_DIR). Not timed: the times of this mode come
# from timing_exclusive.sbatch (arm ref). Ran 2026-11-02 (job 568514) with StrainDesign ~/sd-scip-pricing; set
# SD_TREE=$HOME/sd-kernel-nonzeros for StrainDesign 1.19.2.
#   sbatch kinetic_variants/refmode_modules.sbatch        # from benchmarks/
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
SD=${SD_TREE:-$HOME/sd-scip-pricing}
export PYTHONPATH=$SD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=${REFMOD_OUT:-/scratch/sadreyev_lab/pschneider/kimonu_trip_2026-11-02/refmode_modules}
mkdir -p $OUT
echo "$(hostname) | kimonu $(git -C .. describe --always --dirty) | straindesign $(git -C $SD rev-parse --short HEAD) | $(date -Is)"
ls unblocked/*_fixed.crn.json | xargs -n1 basename | sed 's/_fixed.crn.json//' \
  | xargs -P 6 -I{} $HOME/miniconda3/envs/straindesign/bin/python run_kimonu_refmode.py {} $OUT/kimonu_{}.json
