#!/usr/bin/env bash
#SBATCH --job-name=kimonu_equiv
#SBATCH --partition=sadreyevlab
#SBATCH --exclude=node05,node06,node07,node08
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00
# Equivalence of two Kimonu commits on the 33 networks (refactor_equivalence.py): every network dumped with the src of
# both commits, eight processes at a time with one BLAS thread each, largest networks first, then compared exactly.
# Submit from the scripts directory of a checkout, with the log on scratch:
#   sbatch -o <dir>/%x_%j.log -e <dir>/%x_%j.err kinetic_variants/refactor_equivalence.sbatch <src A> <src B> <out dir>
cd "$SLURM_SUBMIT_DIR"
export SRC_A=$1 SRC_B=$2 OUT=$3
export SD=${SD_TREE:-$HOME/sd-kernel-nonzeros} PY=$HOME/miniconda3/envs/straindesign/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(hostname) | A $(git -C "$SRC_A" rev-parse --short HEAD) $SRC_A | B $(git -C "$SRC_B" rev-parse --short HEAD) $SRC_B" \
     "| straindesign $(git -C "$SD" rev-parse --short HEAD) | $(date -Is)"
mkdir -p "$OUT/A" "$OUT/B"
models=$(for f in $(ls -S "$($PY records.py read UNBLOCKED)"/*_fixed.crn.json); do
           m=$(basename "$f" _fixed.crn.json); grep -qx "$m" kinetic_variants/models33.txt && echo "$m"; done)
for m in $models; do echo "A $m"; echo "B $m"; done |
  xargs -P 8 -L 1 bash -c 'src=SRC_$0; PYTHONPATH=${!src}:$SD $PY -u kinetic_variants/refactor_equivalence.py dump "$1" "$OUT/$0/$1.pkl"' ||
  echo "some dumps failed"
echo "dumps done $(date -Is)"
$PY kinetic_variants/refactor_equivalence.py compare "$OUT/A" "$OUT/B"
echo "exit $? $(date -Is)"
