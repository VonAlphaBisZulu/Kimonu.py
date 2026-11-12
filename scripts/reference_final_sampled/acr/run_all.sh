# three chains at a time over every sampled stem
cd "$(dirname "$0")"; ls Results/concordant/*.mat | xargs -n1 basename | sed 's/.mat$//' | xargs -P3 -I{} bash run_chain.sh {} > chain_all.log 2>&1; echo ALLDONE >> chain_all.log
