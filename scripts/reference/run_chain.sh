#!/usr/bin/env bash
# The R stages of the reference (kinetic modules, then the ACR and ACRR tables) on its concordance results, with the R
# code of the fork at 8a6ed17 (reference/fetch_upstream.sh), run unmodified:
#   bash reference/run_chain.sh exhaustive <model or stem> ...    # results/reference/exhaustive/acr of KIMONU_OUT
#   bash reference/run_chain.sh sampled <model>_seed<S> ...       # results/reference/sampled/acr
# A model name of the exhaustive arm stands for the stem <model>_concordant_fixed. The input is
# acr/Results/concordant/<stem>.mat (ref_finish.sbatch puts it there) with the coupling matrix of
# acr/Results/stoichiometric_coupling/ (taken from networks/reference_coupling_matrices/ of KIMONU_RECORDS when it is
# missing, the same matrix for every stem of a network). The logs <stem>.kin.log and <stem>.acr.log go to acr/, the tables and modules_<stem>.json to acr/Results/.
# The R code expects its results in ../../Results of Code/kinetic_modules, so it runs from a temporary folder whose
# Code/kinetic_modules links the fork's two R files and whose Results links acr/Results.
set -u; export PATH=$HOME/miniconda3/envs/R/bin:$PATH
H=$(cd "$(dirname "$0")" && pwd); PY=$HOME/miniconda3/envs/straindesign/bin/python
U=${UPSTREAM_ALGORITHM:-$H/Upstream_Algorithm}
test -f "$U/Code/kinetic_modules/code_kineticModule_analysis.R" || { echo "no reference code in $U: bash $H/fetch_upstream.sh"; exit 1; }
mode=${1:?exhaustive or sampled}; shift
case $mode in exhaustive) W=$($PY $H/../records.py write REF_EXHAUSTIVE acr) ;; sampled) W=$($PY $H/../records.py write REF_SAMPLED acr) ;;
  *) echo "mode is exhaustive or sampled, not $mode"; exit 1 ;; esac
CM=$($PY $H/../records.py read COUPLING_MATRICES)
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
mkdir -p "$T/Code/kinetic_modules" "$W/Results/concordant" "$W/Results/stoichiometric_coupling" "$W/Results/MetSingle" \
         "$W/Results/MetDouble" "$W/Results/Overview" "$W/Results/Reactions_Giant"
ln -s "$U"/Code/kinetic_modules/*.R "$T/Code/kinetic_modules/"; ln -s "$W/Results" "$T/Results"
cd "$T/Code/kinetic_modules"
for m in "$@"; do case $mode:$m in exhaustive:*_concordant_fixed*|sampled:*) s=$m ;; *) s=${m}_concordant_fixed ;; esac
  c=$(echo "$s" | sed 's/concordant/stoichiometric_coupling/'); net=$(echo "$s" | sed -E 's/_(concordant|seed).*//')
  [ -e "$W/Results/stoichiometric_coupling/$c.mat" ] || cp -p "$CM/${net}_stoichiometric_coupling_fixed.mat" "$W/Results/stoichiometric_coupling/$c.mat"
  echo "== kinetic $m $(date -Is)"; timeout 7200 Rscript code_kineticModule_analysis.R $s.mat > $W/$m.kin.log 2>&1 && echo ok || echo FAIL
  echo "== tables  $m $(date -Is)"; timeout 7200 Rscript Get_excel_tables_kinetic_modules.R $s.RData > $W/$m.acr.log 2>&1 && echo ok || echo FAIL
  Rscript -e "load('$W/Results/$s.RData'); library(jsonlite); writeLines(toJSON(list(modules=lapply(final_lres, function(x) as.integer(x)-1L), balanced=as.integer(s\$B)-1L)), '$W/Results/modules_$s.json')" > /dev/null 2>&1 && echo "modules dumped" || echo "dump FAIL"
done; echo ALLDONE
