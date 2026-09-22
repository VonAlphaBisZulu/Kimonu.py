# Reference R stages (fixed kinetic code from the fork) on the final exhaustive .mat files.
# Run from this directory: bash run_chain.sh <model> [<model> ...]
set -u; export PATH=$HOME/miniconda3/envs/R/bin:$PATH; W=$(cd "$(dirname "$0")" && pwd)
cd $W/Code/kinetic_modules
for m in "$@"; do case $m in *_concordant_fixed*) s=$m;; *) s=${m}_concordant_fixed;; esac
  echo "== kinetic $m $(date -Is)"; timeout 7200 Rscript code_kineticModule_analysis.R $s.mat > $W/$m.kin.log 2>&1 && echo ok || echo FAIL
  echo "== tables  $m $(date -Is)"; timeout 7200 Rscript Get_excel_tables_kinetic_modules.R $s.RData > $W/$m.acr.log 2>&1 && echo ok || echo FAIL
  Rscript -e "load('$W/Results/$s.RData'); library(jsonlite); writeLines(toJSON(list(modules=lapply(final_lres, function(x) as.integer(x)-1L), balanced=as.integer(s\$B)-1L)), '$W/Results/modules_$s.json')" > /dev/null 2>&1 && echo "modules dumped" || echo "dump FAIL"
done; echo ALLDONE
