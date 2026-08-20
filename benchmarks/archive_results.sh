#!/usr/bin/env bash
# Archive the full per-run outputs, not just the summary numbers.
#
# The summary tables are derived; the sets are not. A disagreement in ACRR or in a module
# partition can only be re-examined later if the actual sets were kept, and reruns are expensive
# (COCOA is hours per model and is not reproducible at a fixed seed, so a discarded run cannot be
# recreated even in principle).
set -u
BASE=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks
DEST=${DEST:-$BASE/archive}
STAMP=$(date +%Y%m%d_%H%M)
mkdir -p $DEST
TAR=$DEST/results_${STAMP}.tar.zst

cd $BASE
# Collect only the paths that exist; an empty glob otherwise makes tar exit nonzero and the
# fallback fire on top of a perfectly good archive.
FILES=()
for pat in triplicate_final/*.json sweep/cocoa/cocoa_*.json reference/*.json \
           unblocked/*.crn.json BENCHMARK_DESIGN.md; do
  for f in $pat; do [ -e "$f" ] && FILES+=("$f"); done
done
[ ${#FILES[@]} -eq 0 ] && { echo "nothing to archive"; exit 1; }

if command -v zstd >/dev/null 2>&1; then
  TAR=$DEST/results_${STAMP}.tar.zst
  tar --use-compress-program="zstd -19 -T4" -cf "$TAR" "${FILES[@]}"
else
  TAR=$DEST/results_${STAMP}.tar.gz
  tar -czf "$TAR" "${FILES[@]}"
fi

ls -la "$TAR"
echo "${#FILES[@]} files archived:"
printf '%s\n' "${FILES[@]}" | sed 's|/.*||' | sort | uniq -c
