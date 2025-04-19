library(dplyr)
source("dgps/sim_data_3.R")
truth <- get_truth()

res_df <- read.csv("out/dgp_3_2_0419_141313.csv")

res_df %>%
  summarize(mse_relax = mean((psi_relax - truth)^2),
            mse_tmle = mean((psi_tmle - truth)^2),
            coverage_relax = mean(truth >= lower_relax & truth <= upper_relax),
            coverage_tmle = mean(truth >= lower_tmle & truth <= upper_tmle))
