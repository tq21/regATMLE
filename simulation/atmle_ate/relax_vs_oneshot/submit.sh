#!/bin/bash
# Job name:
#SBATCH --job-name=relax_vs_tmle
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=72:00:00
#
# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=10
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=sky.qiu@berkeley.edu

R CMD BATCH --no-save run_gamma_0.5.R logs/run_gamma_0.5.Rout &
R CMD BATCH --no-save run_gamma_3.R logs/run_gamma_3.Rout &
R CMD BATCH --no-save run_gamma_5.R logs/run_gamma_5.Rout &

wait
