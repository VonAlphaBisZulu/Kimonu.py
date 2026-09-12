#!/usr/bin/env bash
#SBATCH --job-name=cocoa110
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/sweep/cocoa_eff/v110_%A_%a.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
cd $BASE
OUT=benchmarks/sweep/cocoa_eff
mapfile -t CELLS < $BASE/benchmarks/eff_cells.txt
spec=${CELLS[$SLURM_ARRAY_TASK_ID]}
m=${spec%%:*}; SEED=${spec##*:}
dst=$OUT/cocoa_eff_${m}_HiGHS_p64_seed${SEED}.json
[ -s "$dst" ] && { echo "have $m seed $SEED"; exit 0; }
echo "=== $m seed $SEED | COCOA v1.1.0 | HiGHS | 64 workers | $(date -Is)"
/usr/bin/time -v env COCOA_SOLVER=HiGHS COCOA_PROCS=64 COCOA_EFFICIENT=true \
  julia --project=benchmarks/julia_v110 benchmarks/run_cocoa_acr.jl \
    benchmarks/unblocked/${m}_fixed.crn.json $dst $SEED \
  > $OUT/cell_${m}_seed${SEED}.log 2>&1 || echo "cell failed"
grep -E "Elapsed \(wall|Maximum resident|ERROR" $OUT/cell_${m}_seed${SEED}.log || true
echo "--- $m seed $SEED done $(date -Is) ---"
