#!/usr/bin/env bash
#SBATCH --job-name=kimonu_time
#SBATCH --partition=sadreyevlab
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_timing/logs/%x_%j.log
# Runtime benchmark, nothing else on the node: one process at a time, pinned to one core, one BLAS thread.
# Arms ($@): det (iterated route, detached), att (attached), ref (reference-equivalent mode), tan and tanatt (tangent
# route, detached and attached), small (tangent route on the small networks and the ladder). Three runs per network, with
# seeds 0, 1, 2 of the interior flux r*; the reference-equivalent mode has no r* and repeats one computation. A cell whose
# output exists is skipped, so a rerun picks up where a job stopped. Each arm's directory gets the job's provenance and
# the node's hardware (kinetic_variants/provenance.sh). The nodes of partition sadreyevlab are Xeon Silver 4314.
#   sbatch [--mem=2G] timing_exclusive.sbatch det att ref
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
SD=${SD_TREE:-$HOME/sd-kernel-nonzeros}
export PYTHONPATH=$PWD:$SD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=${TIMING_OUT:-/scratch/sadreyev_lab/pschneider/kimonu_timing}
PY=$HOME/miniconda3/envs/straindesign/bin/python
CORE=$($PY -c "import os; print(min(os.sched_getaffinity(0)))")
echo "$(hostname) core $CORE of $(nproc) | kimonu $(git -C .. rev-parse --short HEAD) | straindesign $(git -C $SD rev-parse --short HEAD) | $(date -Is)"
squeue -h -w "$(hostname)" -o "%i %u %j" | head -5
for arm in "$@"; do bash kinetic_variants/provenance.sh "$OUT/$arm" hw; done
models() { for c in unblocked/*_fixed.crn.json; do basename "$c" _fixed.crn.json; done; }
cell() {  # dest, command...
  local dest=$1; shift
  [ -s "$dest" ] && return
  mkdir -p "$(dirname "$dest")"
  /usr/bin/time -v taskset -c "$CORE" "$@" > "${dest%.json}.log" 2>&1
  echo "$(date +%T) rc=$? $(grep -E 'Elapsed|Maximum res' "${dest%.json}.log" | awk '{print $NF}' | tr '\n' ' ') $(basename "$dest")"
}
for R in 1 2 3; do
  SEED=$((R - 1))
  for arm in "$@"; do
    case $arm in
      det|att)
        [ "$arm" = det ] && K=--detach-uptakes || K=
        for m in $(models); do cell $OUT/$arm/kimonu_${m}_rep$R.json $PY run_kimonu_crn.py unblocked/${m}_fixed.crn.json $OUT/$arm/kimonu_${m}_rep$R.json --repeats 1 --assume-unblocked --solver scip --seed $SEED $K; done ;;
      ref)
        for m in $(models); do cell $OUT/$arm/kimonu_${m}_rep$R.json $PY run_kimonu_refmode.py "$m" $OUT/$arm/kimonu_${m}_rep$R.json; done ;;
      tan)
        for m in $(models); do cell $OUT/$arm/tangent_${m}_rep$R.json $PY kinetic_variants/tangent_cell.py "$m" $OUT/$arm/tangent_${m}_rep$R.json --seed $SEED; done ;;
      tanatt)
        for m in $(models); do cell $OUT/$arm/tangent_${m}_rep$R.json $PY kinetic_variants/tangent_cell.py "$m" $OUT/$arm/tangent_${m}_rep$R.json --attached --seed $SEED; done ;;
      small)
        for n in $($PY -c "
import sys; sys.argv=['x']
from kinetic_variants import local_acr_symbolic as L
print(*L.NETWORKS)" 2>/dev/null | tail -1) $(ls /scratch/sadreyev_lab/pschneider/kimonu_local_acr/ladder_*.json | xargs -n1 basename | sed 's/ladder_//; s/.json//'); do
          cell $OUT/$arm/tangent_${n}_rep$R.json $PY kinetic_variants/tangent_cell.py "$n" $OUT/$arm/tangent_${n}_rep$R.json --attached --seed $SEED
        done ;;
    esac
  done
done
echo "done $(date -Is)"
