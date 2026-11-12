#!/usr/bin/env bash
#SBATCH --job-name=refmode-mod
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-02:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_trip_2026-11-02/logs/refmod_%j.log
# Kinetic modules of the reference-equivalent mode on all 33 networks, six networks at a time, for the coupled pairs
# of Figure 4C, D in early November, written to experiments/refmode_modules/ of KIMONU_OUT (records.py). Superseded: pair_sets.py
# reads the modules of the timing arm ref (timing_exclusive.sbatch), whose runs give the times of this mode too. Ran 2026-11-02 (job 568514) with StrainDesign ~/sd-scip-pricing; set
# SD_TREE=$HOME/sd-kernel-nonzeros for StrainDesign 1.19.2.
#   sbatch kinetic_variants/refmode_modules.sbatch        # from scripts/
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts
SD=${SD_TREE:-$HOME/sd-scip-pricing}
export PYTHONPATH=$SD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=$($HOME/miniconda3/envs/straindesign/bin/python records.py write EXPERIMENTS refmode_modules)
echo "$(hostname) | kimonu $(git -C .. describe --always --dirty) | straindesign $(git -C $SD rev-parse --short HEAD) | $(date -Is)"
ls "$($HOME/miniconda3/envs/straindesign/bin/python records.py read UNBLOCKED)"/*_fixed.crn.json | xargs -n1 basename | sed 's/_fixed.crn.json//' \
  | xargs -P 6 -I{} $HOME/miniconda3/envs/straindesign/bin/python run_kimonu_refmode.py {} $OUT/kimonu_{}.json
