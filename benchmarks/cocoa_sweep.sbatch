#!/usr/bin/env bash
#SBATCH --job-name=cocoasweep
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/sweep/cocoa/sweep_%j.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
cd $BASE
OUT=benchmarks/sweep/cocoa
SEED=${SEED:-1234}
SKIP=${SKIP:-}

# Smallest first: COCOA's cost grows with complexes^2, so this banks results early rather than
# stalling on a model that may not finish inside the walltime.
for f in $(ls -S -r benchmarks/unblocked/*_fixed.crn.json); do
  m=$(basename $f _fixed.crn.json)
  case " $SKIP " in *" $m "*) echo "skip $m"; continue;; esac
  dst=$OUT/cocoa_unb_${m}_HiGHS_p64_seed${SEED}.json
  [ -s "$dst" ] && { echo "have $m already"; continue; }
  echo "=== $m | HiGHS | 64 workers | $(date -Is)"
  CELL=$OUT/cell_unb_${m}.log
  /usr/bin/time -v env COCOA_SOLVER=HiGHS COCOA_PROCS=64 \
    julia --project=benchmarks/julia_revision benchmarks/run_cocoa_acr.jl \
      benchmarks/unblocked/${m}_fixed.crn.json $dst $SEED \
    > $CELL 2>&1 || echo "cell failed, see $CELL"
  grep -E "Elapsed time|seconds =|Maximum resident|Elapsed \(wall|ERROR" $CELL || true
  echo "--- $m done $(date -Is) ---"
done
echo "ALL COCOA CELLS DONE $(date -Is)"
