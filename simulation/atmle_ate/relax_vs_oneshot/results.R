library(dplyr)
source("dgps/sim_data_2.R")
truth <- get_truth()

res_df <- read.csv("out/dgp_2_0.5_0608_130212.csv")

# TODO: add lambda

res_df %>%
  filter(j == 5) %>%
  summarize(abs_bias_relax = abs(mean(psi_relax - truth)),
            abs_bias_tmle = abs(mean(psi_tmle - truth)),
            se_relax = sd(psi_relax),
            se_tmle = sd(psi_tmle),
            mse_relax = mean((psi_relax - truth)^2),
            mse_tmle = mean((psi_tmle - truth)^2),
            cover_relax = mean(truth >= lower_relax & truth <= upper_relax),
            cover_tmle = mean(truth >= lower_tmle & truth <= upper_tmle),
            oracle_cover_relax = mean(truth >= psi_relax - 1.96 * sd(psi_relax) & truth <= psi_relax + 1.96 * sd(psi_relax)),
            oracle_cover_tmle = mean(truth >= psi_tmle - 1.96 * sd(psi_tmle) & truth <= psi_tmle + 1.96 * sd(psi_tmle)),
            .by = "n")
