#!/usr/bin/env bash
#SBATCH --job-name=cocoa5k
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts/sweep/cocoa_s5000/s5k_%A_%a.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
cd $BASE
OUT=scripts/sweep/cocoa_s5000
mapfile -t CELLS < $BASE/scripts/s5000_cells.txt
spec=${CELLS[$SLURM_ARRAY_TASK_ID]}
m=${spec%%:*}; SEED=${spec##*:}
dst=$OUT/cocoa_s5000_${m}_HiGHS_p64_seed${SEED}.json
[ -s "$dst" ] && { echo "have $m seed $SEED"; exit 0; }
echo "=== $m seed $SEED | v1.1.0 | sample_size=5000 | $(date -Is)"
/usr/bin/time -v env COCOA_SOLVER=HiGHS COCOA_PROCS=64 COCOA_SAMPLES=5000 \
  julia --project=scripts/julia_v110 scripts/run_cocoa_acr.jl \
    scripts/unblocked/${m}_fixed.crn.json $dst $SEED \
  > $OUT/cell_${m}_seed${SEED}.log 2>&1 || echo "cell failed"
grep -E "Elapsed \(wall|ERROR" $OUT/cell_${m}_seed${SEED}.log || true
echo "--- $m seed $SEED done $(date -Is) ---"
