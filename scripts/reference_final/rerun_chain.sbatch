#!/usr/bin/env bash
#SBATCH --job-name=ref-chain
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts/reference_final/chain_%j.log
# R stages again for stems whose coupling matrix was missing: bash run_chain.sh <stem> ...
bash /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts/reference_final/acr/run_chain.sh "$@"
