#!/usr/bin/env bash
# Rerun the reference's balanced-complex stage on its own published models and compare against the
# committed results, while recording the sign-class census. The census bounds how much of each
# network depends on the negative branch: any complex classified N or U can only be tested as a
# denominator through that branch, whose LP needed adjusting at genome scale (pull request 1 of the fork).
set -u
UPSTREAM=${UPSTREAM_ALGORITHM:?set UPSTREAM_ALGORITHM}
HERE=$(cd "$(dirname "$0")" && pwd)
for name in "$@"; do
  model="$UPSTREAM/Models/models_with_elementary_steps/${name}_pre_balanced_fixed.mat"
  published="$UPSTREAM/Results/concordant/${name}_concordant_fixed.mat"
  if [ ! -f "$model" ] || [ ! -f "$published" ]; then
    echo "== $name: missing model or published result, skipping"
    continue
  fi
  echo "== $name"
  octave --no-gui "$HERE/validate_against_published.m" "$model" "$published" 2>&1 \
    | grep -E "complexes,|rerun:|sign classes|exposed|published:|identical:"
done
