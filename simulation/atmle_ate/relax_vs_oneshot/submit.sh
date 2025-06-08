#!/bin/bash
# Job name:
#SBATCH --job-name=relax_vs_tmle
#
# Partition:
#SBATCH --partition=savio3
#
#SBATCH --qos=biostat_savio3_normal
#SBATCH --account=co_biostat
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=72:00:00
#
# Number of nodes for use case:
#SBATCH --nodes=7
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=10
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=sky.qiu@berkeley.edu

module load r

R CMD BATCH --no-save run_dgp_1.R logs/run_dgp_1.Rout &
R CMD BATCH --no-save run_dgp_2_gamma_0.5.R logs/run_dgp_2_gamma_0.5.Rout &
R CMD BATCH --no-save run_dgp_2_gamma_1.R logs/run_dgp_2_gamma_1.Rout &
R CMD BATCH --no-save run_dgp_2_gamma_2.R logs/run_dgp_2_gamma_2.Rout &
R CMD BATCH --no-save run_dgp_3_gamma_0.5.R logs/run_dgp_3_gamma_0.5.Rout &
R CMD BATCH --no-save run_dgp_3_gamma_1.R logs/run_dgp_3_gamma_1.Rout &
R CMD BATCH --no-save run_dgp_3_gamma_2.R logs/run_dgp_3_gamma_2.Rout &

wait
