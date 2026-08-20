#!/usr/bin/env bash
#SBATCH --job-name=cocoadet
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/sweep/cocoa/determinism_%j.log
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
OUT=benchmarks/sweep/cocoa/determinism
mkdir -p $OUT
M=${M:-iCN900}

# Same model, same seed, twice at each worker count. If the serial pair agrees and the parallel
# pair does not, the non-determinism is in how pmap splits candidate pairs across workers, not in
# the RNG. iCN900 is used rather than iLJ478 because it is ~15x cheaper and showed the same shape.
for procs in 1 1 64 64; do
  tag=$OUT/${M}_p${procs}_$(date +%s%N).json
  echo "=== $M procs=$procs seed=1234 $(date -Is)"
  env COCOA_SOLVER=HiGHS COCOA_PROCS=$procs \
    julia --project=benchmarks/julia_revision benchmarks/run_cocoa_acr.jl \
      benchmarks/unblocked/${M}_fixed.crn.json $tag 1234 > ${tag%.json}.log 2>&1 \
    || echo "failed, see ${tag%.json}.log"
  python3 -c "
import json,sys
try:
    d=json.load(open('$tag'))
    print(f'  procs=$procs -> ACR {len(d[\"acr\"])} ACRR {len(d[\"acrr\"])} giant {d[\"giant_module\"]} secs {d[\"seconds\"]:.0f}')
except Exception as e: print('  no output:', e)"
done
echo "DONE $(date -Is)"
