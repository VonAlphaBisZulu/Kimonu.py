#!/usr/bin/env bash
# Records the code and software of a run in its output directory, as provenance_<job>.txt, and with "hw" also the
# node's hardware (hostname, lscpu, free -g) as hardware_<job>.txt. An array writes once, from its first task.
#   bash kinetic_variants/provenance.sh <output directory> [hw]
D=$1; mkdir -p "$D"
[ -n "$SLURM_ARRAY_TASK_ID" ] && [ "$SLURM_ARRAY_TASK_ID" != "${SLURM_ARRAY_TASK_MIN:-$SLURM_ARRAY_TASK_ID}" ] && exit 0
J=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}
K=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SD=${SD_TREE:-$HOME/sd-kernel-nonzeros}
PY=$HOME/miniconda3/envs/straindesign/bin/python
{
  echo "date $(date -Is)"
  echo "host $(hostname)"
  echo "job $J"
  echo "kimonu $(git -C "$K" rev-parse HEAD)"
  echo "kimonu_uncommitted_files $(git -C "$K" status --porcelain --untracked-files=no -- src scripts | wc -l)"
  echo "straindesign_tree $SD"
  echo "straindesign $(git -C "$SD" describe --tags --always) $(git -C "$SD" rev-parse HEAD)"
  echo "python $($PY -V 2>&1)"
} > "$D/provenance_$J.txt"
[ "$2" = hw ] && { hostname; echo; lscpu; echo; free -g; } > "$D/hardware_$J.txt"
exit 0
