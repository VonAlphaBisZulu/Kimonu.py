#!/usr/bin/env bash
#SBATCH --job-name=cocoa110
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_records/results/cocoa_jl/efficient_true/v110_%A_%a.log
# Frozen harness of the paper (COCOA.jl v1.1.0). Outputs go to results/cocoa_jl/efficient_true/ of KIMONU_OUT (records.py); the log line
# above names the records folder of the cluster, so with another KIMONU_OUT submit with -o <that folder>/v110_%A_%a.log.
set -u
export JULIA_DEPOT_PATH=/scratch/sadreyev_lab/pschneider/.julia
export PATH=$HOME/.local/bin:$PATH
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
cd $BASE
PY=$HOME/miniconda3/envs/straindesign/bin/python
OUT=$($PY scripts/records.py write COCOA_TRUE)
NET=$($PY scripts/records.py read UNBLOCKED)
mapfile -t CELLS < $BASE/scripts/eff_big_cells.txt
spec=${CELLS[$SLURM_ARRAY_TASK_ID]}
m=${spec%%:*}; SEED=${spec##*:}
dst=$OUT/cocoa_eff_${m}_HiGHS_p64_seed${SEED}.json
[ -s "$dst" ] && { echo "have $m seed $SEED"; exit 0; }
echo "=== $m seed $SEED | COCOA v1.1.0 | HiGHS | 64 workers | $(date -Is)"
/usr/bin/time -v env COCOA_SOLVER=HiGHS COCOA_PROCS=64 COCOA_EFFICIENT=true \
  julia --project=scripts/julia_v110 scripts/run_cocoa_acr.jl \
    $NET/${m}_fixed.crn.json $dst $SEED \
  > $OUT/cell_${m}_seed${SEED}.log 2>&1 || echo "cell failed"
grep -E "Elapsed \(wall|Maximum resident|ERROR" $OUT/cell_${m}_seed${SEED}.log || true
echo "--- $m seed $SEED done $(date -Is) ---"
