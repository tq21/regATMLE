library(dplyr)
source("sim_data.R")
truth <- get_truth()

res_df <- read.csv("out/res_0312_133139.csv")

tmp <- res_df %>%
  summarize(mse_proj = mean((psi_proj - truth)^2),
            #mse_proj_cv = mean((psi_proj_cv - truth)^2),
            mse_delta = mean((psi_delta - truth)^2),
            coverage_proj = mean(truth >= lower_proj & truth <= upper_proj),
            #coverage_proj_cv = mean(truth >= lower_proj_cv & truth <= upper_proj_cv),
            coverage_delta = mean(truth >= lower_delta & truth <= upper_delta),
            .by = "j")

sum(res_df$j == 2)
