#!/usr/bin/env bash
# Fetches the reference implementation of Langary et al. (2025), Upstream_Algorithm, at the commit the paper ran: the
# head 8a6ed17 of the branch fix/negative-branch-charnes-cooper of the fork https://github.com/VonAlphaBisZulu/Upstream_Algorithm,
# which is pull request 1 of https://github.com/ankueken/Upstream_Algorithm (the settings that needed adjusting at
# genome scale, Supplementary Text S3). The clone is not part of this repository and is ignored by git.
#   bash reference/fetch_upstream.sh [--full] [<commit> [<dir>]]      # default 8a6ed17 into scripts/reference/Upstream_Algorithm
# By default only the code (Code/ with the MATLAB and R files, and the top-level files) is checked out, from a partial
# clone of a few MB, which is all that run_reference.m, retest_rejected.m and reference/run_chain.sh need. --full checks
# out the models and published results as well (Models/, Results/, about 2 GB), which export_split_model.py,
# compare_to_published.py, validate_*.py, matlab_reference/census.sh and kinetic_variants/lumped_toric.py read; put
# that clone on scratch. The scripts that run the reference take the clone from UPSTREAM_ALGORITHM, default the folder
# above. Another commit, for example the original code of ankueken/Upstream_Algorithm at d729161, goes to its own folder:
#   bash reference/fetch_upstream.sh d729161 /tmp/Upstream_Algorithm_d729161
set -eu
H=$(cd "$(dirname "$0")" && pwd)
FULL=0; [ "${1:-}" = --full ] && { FULL=1; shift; }
PR_HEAD=8a6ed17384812f1a5b973a0752174b6a7797e8fa
C=${1:-$PR_HEAD}; D=${2:-$H/Upstream_Algorithm}
URL=https://github.com/VonAlphaBisZulu/Upstream_Algorithm.git
if [ ! -d "$D/.git" ]; then
  git clone --quiet --filter=blob:none --no-checkout "$URL" "$D"
  git -C "$D" remote add upstream https://github.com/ankueken/Upstream_Algorithm.git
fi
git -C "$D" fetch --quiet --filter=blob:none origin fix/negative-branch-charnes-cooper
git -C "$D" fetch --quiet --filter=blob:none upstream main 2>/dev/null || true
if [ $FULL = 1 ]; then git -C "$D" sparse-checkout disable; else git -C "$D" sparse-checkout set --cone Code; fi
git -C "$D" -c advice.detachedHead=false checkout --quiet "$C"
got=$(git -C "$D" rev-parse HEAD)
case $got in "$C"*) ;; *) echo "checked out $got, not $C" >&2; exit 1 ;; esac
[ -z "$(git -C "$D" status --porcelain --untracked-files=no)" ] || { echo "$D has local changes" >&2; exit 1; }
echo "Upstream_Algorithm at $(git -C "$D" log -1 --format='%h %ad %s' --date=short) in $D ($([ $FULL = 1 ] && echo full || echo code only))"
