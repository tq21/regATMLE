#!/bin/bash
# Job name:
#SBATCH --job-name=atmle_cv_wm
#
# Partition:
#SBATCH --partition=savio3
#
#SBATCH --qos=biostat_savio3_normal
#SBATCH --account=co_biostat
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=10:00:00
#
# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=yi_li@berkeley.edu

module load r

R CMD BATCH --no-save condition_number_comparison.R logs/condition_number_comparison.Rout &

wait
