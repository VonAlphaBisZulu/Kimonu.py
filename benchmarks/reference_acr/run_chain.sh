set -u; W=/tmp/claude-1031/acr; export PATH=$HOME/miniconda3/envs/R/bin:$PATH
cd $W/in_patched
for f in *_concordant_fixed.mat; do echo "== kinetic $f $(date -Is)"; timeout 3600 Rscript km_fixed.R $f > ../${f%.mat}.kin.log 2>&1 && echo "ok kinetic $f" || echo "FAIL kinetic $f"; done
# km_fixed.R saves ../<name>.RData i.e. into $W ; move into Results/ where the ACR script loads from
mv $W/*_concordant_fixed.RData $W/Results/ 2>/dev/null
cd $W/Code/kinetic_modules
for f in ../../Results/*_concordant_fixed.RData; do kf=$(basename $f); echo "== acr $kf $(date -Is)"; timeout 3600 Rscript Get_excel_tables_kinetic_modules.R $kf > ../../${kf%.RData}.acr.log 2>&1 && echo "ok acr $kf" || echo "FAIL acr $kf"; done
echo ALLDONE
