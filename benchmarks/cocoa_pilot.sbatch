#!/usr/bin/env bash
#SBATCH --job-name=cocoapilot
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/sweep/cocoa/pilot_%j.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=benchmarks/sweep/cocoa
for m in iIS312_Trypomastigote iIS312_Amastigote iLJ478; do
  echo "=== $m | HiGHS | 64 workers | $(date -Is)"
  CELL=$OUT/cell_pilot_${m}.log
  /usr/bin/time -v env COCOA_SOLVER=HiGHS COCOA_PROCS=64 \
    julia --project=benchmarks/julia_revision benchmarks/run_cocoa_acr.jl \
      $OUT/${m}_fixed.crn.json $OUT/cocoa_${m}_HiGHS_p64_seed1234.json 1234 \
    > $CELL 2>&1 || echo "cell failed, see $CELL"
  grep -E "Elapsed time|seconds =|Maximum resident|Elapsed \(wall|ERROR" $CELL || true
  echo "--- cell done: $(date -Is) ---"
done
