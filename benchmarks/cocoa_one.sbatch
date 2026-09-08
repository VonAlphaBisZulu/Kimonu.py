#!/usr/bin/env bash
#SBATCH --job-name=cocoa1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/sweep/cocoa/one_%j.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
cd $BASE
OUT=benchmarks/sweep/cocoa
SEED=${SEED:-1234}
m=$MODEL
dst=$OUT/cocoa_unb_${m}_HiGHS_p64_seed${SEED}.json
[ -s "$dst" ] && { echo "have $m already"; exit 0; }
echo "=== $m | HiGHS | 64 workers | $(date -Is)"
# Same per-model configuration as the sweep -- full node, 64 workers -- so wall times stay
# comparable across models; only the assignment of models to nodes changes.
/usr/bin/time -v env COCOA_SOLVER=HiGHS COCOA_PROCS=64 \
  julia --project=benchmarks/julia_revision benchmarks/run_cocoa_acr.jl \
    benchmarks/unblocked/${m}_fixed.crn.json $dst $SEED \
  > $OUT/cell_unb_${m}.log 2>&1 || echo "cell failed, see $OUT/cell_unb_${m}.log"
grep -E "Elapsed \(wall|Maximum resident|ERROR" $OUT/cell_unb_${m}.log || true
echo "--- $m done $(date -Is) ---"
