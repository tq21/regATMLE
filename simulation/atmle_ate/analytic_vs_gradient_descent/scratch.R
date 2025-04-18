library(dplyr)
truth <- get_truth()

res_df %>%
  summarize(mse_relax_an = mean((psi_relax_an - truth)^2),
            mse_relax_gd = mean((psi_relax_gd - truth)^2),
            coverage_relax_an = mean(truth >= lower_relax_an & truth <= upper_relax_an),
            coverage_relax_gd = mean(truth >= lower_relax_gd & truth <= upper_relax_gd))
