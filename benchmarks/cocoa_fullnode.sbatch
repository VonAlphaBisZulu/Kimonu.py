#!/usr/bin/env bash
#SBATCH --job-name=cocoa64
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/sweep/cocoa/fullnode_%j.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
export CPLEX_STUDIO_BINARIES=$HOME/miniconda3/envs/straindesign/lib/python3.12/site-packages/cplex/_internal
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py

MODEL=${MODEL:?set MODEL}
SEEDS=${SEEDS:-"1234 42 7 999 2026 31337"}
PROCS=${PROCS:-64}
SOLVER=${SOLVER:-Gurobi}
OUT=benchmarks/sweep/cocoa

for s in $SEEDS; do
  echo "=== $MODEL seed $s | $SOLVER | $PROCS workers | $(date -Is)"
  # Full output is kept: a filtered log hides licence failures and crashes, and leaves nothing to
  # tell "still running" from "died". The summary is extracted afterwards, not instead.
  CELL=$OUT/cell_${MODEL}_${SOLVER}_p${PROCS}_seed${s}.log
  /usr/bin/time -v env COCOA_SOLVER=$SOLVER COCOA_PROCS=$PROCS \
    julia --project=benchmarks/julia_revision benchmarks/run_cocoa_acr.jl \
      $OUT/${MODEL}_fixed.crn.json $OUT/cocoa_${MODEL}_${SOLVER}_p${PROCS}_seed${s}.json $s \
    > $CELL 2>&1 || echo "cell failed, see $CELL"
  grep -E "Elapsed time|seconds =|Maximum resident|Elapsed \(wall|ERROR" $CELL || true
  echo "--- cell done: $(date -Is) ---"
done
