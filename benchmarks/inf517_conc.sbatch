#!/usr/bin/env bash
#SBATCH --job-name=inf517conc
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=6:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/sweep/modules/inf517conc_%j.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=benchmarks/sweep/modules
mkdir -p $OUT
for M in iNF517 iSB619; do
  echo "=== $M $(date -Is)"
  env COCOA_SOLVER=HiGHS COCOA_PROCS=64 \
    julia --project=benchmarks/julia_revision benchmarks/run_cocoa_acr.jl \
      benchmarks/unblocked/${M}_fixed.crn.json $OUT/cocoa_conc_${M}.json 1234 \
      > $OUT/cocoa_conc_${M}.log 2>&1 || echo "failed $M"
done
echo "DONE $(date -Is)"
