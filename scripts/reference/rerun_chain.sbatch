#!/usr/bin/env bash
#SBATCH --job-name=ref-chain
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_records/results/reference/exhaustive/chain_%j.log
# R stages again for stems whose coupling matrix was missing: sbatch reference/rerun_chain.sbatch <stem> ...
# The log goes to the records folder of the cluster; with another KIMONU_OUT, submit with -o <that folder>/results/reference/exhaustive/chain_%j.log.
bash /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts/reference/run_chain.sh exhaustive "$@"
