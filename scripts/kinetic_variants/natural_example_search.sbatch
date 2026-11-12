#!/usr/bin/env bash
#SBATCH --job-name=natural-toy
#SBATCH --partition=sadreyevlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=0-06:00:00
#SBATCH --output=/scratch/sadreyev_lab/pschneider/kimonu_natural/%x_%j.log
cd /home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts/kinetic_variants
export OMP_NUM_THREADS=8
$HOME/miniconda3/envs/straindesign/bin/python -u natural_example_search.py
