#!/usr/bin/env bash
# The R stages on every sampled stem of results/reference/sampled/acr/Results/concordant/ (KIMONU_RECORDS), three at a
# time, the log in chain_all.log of that folder of KIMONU_OUT:   bash reference/run_all_sampled.sh
H=$(cd "$(dirname "$0")" && pwd); PY=$HOME/miniconda3/envs/straindesign/bin/python
IN=$($PY $H/../records.py read REF_SAMPLED acr Results concordant); W=$($PY $H/../records.py write REF_SAMPLED acr)
ls "$IN"/*.mat | xargs -n1 basename | sed 's/.mat$//' | xargs -P3 -I{} bash "$H/run_chain.sh" sampled {} > "$W/chain_all.log" 2>&1
echo ALLDONE >> "$W/chain_all.log"
