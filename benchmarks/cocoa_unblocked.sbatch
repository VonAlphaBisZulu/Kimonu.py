#!/usr/bin/env bash
#SBATCH --job-name=cocoaunb
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/sweep/cocoa/unblocked_%j.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=benchmarks/sweep/cocoa
echo "=== iJB785 UNBLOCKED | HiGHS | 64 workers | $(date -Is)"
CELL=$OUT/cell_iJB785_unblocked.log
/usr/bin/time -v env COCOA_SOLVER=HiGHS COCOA_PROCS=64 \
  julia --project=benchmarks/julia_revision benchmarks/run_cocoa_acr.jl \
    $OUT/iJB785_unblocked.crn.json $OUT/cocoa_iJB785_unblocked_HiGHS_p64_seed1234.json 1234 \
  > $CELL 2>&1 || echo "cell failed, see $CELL"
grep -E "Elapsed time|seconds =|Maximum resident|Elapsed \(wall|ERROR" $CELL || true
echo "--- done $(date -Is) ---"
