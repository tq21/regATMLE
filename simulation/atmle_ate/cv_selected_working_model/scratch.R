library(dplyr)
source("dgp/sim_data_2.R")
truth <- get_truth_2()

res_df <- read.csv("out/res_dgp_2_0318_175800.csv")

names(res_df)<-c(
  "n", "B", "j",
  "relax_lambda", "proj_lambda", "delta_lambda",
  "psi_relax", "psi_proj", "psi_delta",
  # Non-parametric EIC-based inference
  "se_relax_np", "se_proj_np", "se_delta_np",
  "lower_relax_np", "upper_relax_np",
  "lower_proj_np", "upper_proj_np",
  "lower_delta_np", "upper_delta_np",
  # Projected EIC-based inference(weak)
  "se_relax_proj", "se_proj_proj", "se_delta_proj",
  "lower_relax_proj", "upper_relax_proj",
  "lower_proj_proj", "upper_proj_proj",
  "lower_delta_proj", "upper_delta_proj",
  # Projected EIC-based inference(cv)
  "se_relax_proj_cv", "se_proj_proj_cv", "se_delta_proj_cv",
  "lower_relax_proj_cv", "upper_relax_proj_cv",
  "lower_proj_proj_cv", "upper_proj_proj_cv",
  "lower_delta_proj_cv", "upper_delta_proj_cv",
  # Delta EIC-based inference
  "se_relax_delta", "se_proj_delta", "se_delta_delta",
  "lower_relax_delta", "upper_relax_delta",
  "lower_proj_delta", "upper_proj_delta",
  "lower_delta_delta", "upper_delta_delta"
)

res_df_CLEAN<-na.omit(res_df[,1:45])
res_df_CLEAN %>%filter(j==1) %>%
  summarize(abs_bias_relax = abs(mean(psi_relax - truth)),
            abs_bias_proj = abs(mean(psi_proj - truth)),
            abs_bias_delta = abs(mean(psi_delta - truth)),
            se_relax = sd(psi_relax),
            se_proj = sd(psi_proj),
            se_delta = sd(psi_delta),
            mse_relax = mean((psi_relax - truth)^2),
            mse_proj = mean((psi_proj - truth)^2),
            mse_delta = mean((psi_delta - truth)^2),
            cover_relax_np = mean(truth >= lower_relax_np & truth <= upper_relax_np),
            cover_proj_np = mean(truth >= lower_proj_np & truth <= upper_proj_np),
            cover_delta_np = mean(truth >= lower_delta_np & truth <= upper_delta_np),

            cover_relax_proj = mean(truth >= lower_relax_proj & truth <= upper_relax_proj),
            cover_proj_proj = mean(truth >= lower_proj_proj & truth <= upper_proj_proj),
            cover_delta_proj = mean(truth >= lower_delta_proj & truth <= upper_delta_proj),

            cover_relax_proj_cv = mean(truth >= lower_relax_proj_cv & truth <= upper_relax_proj_cv),
            cover_proj_proj_cv = mean(truth >= lower_proj_proj_cv & truth <= upper_proj_proj_cv),
            cover_delta_proj_cv = mean(truth >= lower_delta_proj_cv & truth <= upper_delta_proj_cv),


            cover_relax_delta = mean(truth >= lower_relax_delta & truth <= upper_relax_delta),
            cover_proj_delta = mean(truth >= lower_proj_delta & truth <= upper_proj_delta),
            cover_delta_delta = mean(truth >= lower_delta_delta & truth <= upper_delta_delta),

            cover_oracle_relax = mean(truth >= psi_relax-1.96*sd(psi_relax) & truth <= psi_relax+1.96*sd(psi_relax)),
            cover_oracle_proj = mean(truth >= psi_proj-1.96*sd(psi_proj) & truth <= psi_proj+1.96*sd(psi_proj)),
            cover_oracle_delta = mean(truth >= psi_delta-1.96*sd(psi_delta) & truth <= psi_delta+1.96*sd(psi_delta)),
            .by = n)


###################################################################################################################

# Load required packages
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

# Compute summary statistics
temp <- res_df_CLEAN %>%
  filter(j == 2) %>%
  group_by(n) %>%
  summarize(
    abs_bias_relax = abs(mean(psi_relax - truth)),
    abs_bias_proj  = abs(mean(psi_proj - truth)),
    abs_bias_delta = abs(mean(psi_delta - truth)),
    se_relax       = sd(psi_relax),
    se_proj        = sd(psi_proj),
    se_delta       = sd(psi_delta),
    mse_relax      = mean((psi_relax - truth)^2),
    mse_proj       = mean((psi_proj - truth)^2),
    mse_delta      = mean((psi_delta - truth)^2),
    # Non-parametric EIC-based inference (coverage and CI length)
    cov_relax_np   = sprintf("%.2f (%.3f)",
                             mean(truth >= lower_relax_np & truth <= upper_relax_np) * 100,
                             mean(upper_relax_np - lower_relax_np)),
    cov_proj_np    = sprintf("%.2f (%.3f)",
                             mean(truth >= lower_proj_np & truth <= upper_proj_np) * 100,
                             mean(upper_proj_np - lower_proj_np)),
    cov_delta_np   = sprintf("%.2f (%.3f)",
                             mean(truth >= lower_delta_np & truth <= upper_delta_np) * 100,
                             mean(upper_delta_np - lower_delta_np)),
    # Projected EIC-based inference (weak)
    cov_relax_proj = sprintf("%.2f (%.3f)",
                             mean(truth >= lower_relax_proj & truth <= upper_relax_proj) * 100,
                             mean(upper_relax_proj - lower_relax_proj)),
    cov_proj_proj  = sprintf("%.2f (%.3f)",
                             mean(truth >= lower_proj_proj & truth <= upper_proj_proj) * 100,
                             mean(upper_proj_proj - lower_proj_proj)),
    cov_delta_proj = sprintf("%.2f (%.3f)",
                             mean(truth >= lower_delta_proj & truth <= upper_delta_proj) * 100,
                             mean(upper_delta_proj - lower_delta_proj)),
    # Projected EIC-based inference (cv)
    cov_relax_proj_cv = sprintf("%.2f (%.3f)",
                                mean(truth >= lower_relax_proj_cv & truth <= upper_relax_proj_cv) * 100,
                                mean(upper_relax_proj_cv - lower_relax_proj_cv)),
    cov_proj_proj_cv  = sprintf("%.2f (%.3f)",
                                mean(truth >= lower_proj_proj_cv & truth <= upper_proj_proj_cv) * 100,
                                mean(upper_proj_proj_cv - lower_proj_proj_cv)),
    cov_delta_proj_cv = sprintf("%.2f (%.3f)",
                                mean(truth >= lower_delta_proj_cv & truth <= upper_delta_proj_cv) * 100,
                                mean(upper_delta_proj_cv - lower_delta_proj_cv)),
    # Delta EIC-based inference
    cov_relax_delta = sprintf("%.2f (%.3f)",
                              mean(truth >= lower_relax_delta & truth <= upper_relax_delta) * 100,
                              mean(upper_relax_delta - lower_relax_delta)),
    cov_proj_delta  = sprintf("%.2f (%.3f)",
                              mean(truth >= lower_proj_delta & truth <= upper_proj_delta) * 100,
                              mean(upper_proj_delta - lower_proj_delta)),
    cov_delta_delta = sprintf("%.2f (%.3f)",
                              mean(truth >= lower_delta_delta & truth <= upper_delta_delta) * 100,
                              mean(upper_delta_delta - lower_delta_delta)),
    # Oracle coverage (CI length computed from 1.96*sd)
    cov_oracle_relax = sprintf("%.2f (%.3f)",
                               mean(truth >= psi_relax - 1.96*sd(psi_relax) & truth <= psi_relax + 1.96*sd(psi_relax)) * 100,
                               2 * 1.96 * sd(psi_relax)),
    cov_oracle_proj  = sprintf("%.2f (%.3f)",
                               mean(truth >= psi_proj - 1.96*sd(psi_proj) & truth <= psi_proj + 1.96*sd(psi_proj)) * 100,
                               2 * 1.96 * sd(psi_proj)),
    cov_oracle_delta = sprintf("%.2f (%.3f)",
                               mean(truth >= psi_delta - 1.96*sd(psi_delta) & truth <= psi_delta + 1.96*sd(psi_delta)) * 100,
                               2 * 1.96 * sd(psi_delta))
  )

# Reshape data for each targeting method
temp_relax <- temp %>% transmute(
  n,
  Targeting = "Relaxed",
  Abs_Bias = abs_bias_relax,
  Std_Err = se_relax,
  MSE = mse_relax,
  `Cov NP (%)` = cov_relax_np,
  `Cov Proj (%)` = cov_relax_proj,
  `Cov Proj CV (%)` = cov_relax_proj_cv,
  `Cov Delta (%)` = cov_relax_delta,
  `Oracle Cov (%)` = cov_oracle_relax
)

temp_proj <- temp %>% transmute(
  n,
  Targeting = "Projection",
  Abs_Bias = abs_bias_proj,
  Std_Err = se_proj,
  MSE = mse_proj,
  `Cov NP (%)` = cov_proj_np,
  `Cov Proj (%)` = cov_proj_proj,
  `Cov Proj CV (%)` = cov_proj_proj_cv,
  `Cov Delta (%)` = cov_proj_delta,
  `Oracle Cov (%)` = cov_oracle_proj
)

temp_delta <- temp %>% transmute(
  n,
  Targeting = "Delta-method",
  Abs_Bias = abs_bias_delta,
  Std_Err = se_delta,
  MSE = mse_delta,
  `Cov NP (%)` = cov_delta_np,
  `Cov Proj (%)` = cov_delta_proj,
  `Cov Proj CV (%)` = cov_delta_proj_cv,
  `Cov Delta (%)` = cov_delta_delta,
  `Oracle Cov (%)` = cov_oracle_delta
)

# Combine the data frames and order by sample size and targeting method
temp_table <- bind_rows(temp_relax, temp_proj, temp_delta) %>%
  arrange(n, factor(Targeting, levels = c("Relaxed", "Projection", "Delta-method")))

# Format numeric columns: Bias and Std. Err. to 4 decimals, MSE to 6 decimals.
temp_table <- temp_table %>%
  mutate(
    Abs_Bias = sprintf("%.4f", Abs_Bias),
    Std_Err  = sprintf("%.4f", Std_Err),
    MSE      = sprintf("%.6f", MSE)
  )

# Add an Estimator column so that "A-TMLE" appears in every row
temp_table <- temp_table %>%
  mutate(Estimator = "A-TMLE") %>%
  select(Estimator, everything())

# Generate the LaTeX table using kable, collapse the Estimator column,
# and add a horizontal line after every 3 rows (from columns 2 to 11)
latex_table <- kable(temp_table, format = "latex", booktabs = TRUE,
                     caption = "Simulation results for A-TMLE. DGP 1. Undersmoothed",
                     align = c("c", "c", "c", "r", "r", "r", "r", "r", "r", "r", "r")) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
  row_spec(3, extra_latex_after = "\\cmidrule(lr){2-11}") %>%
  row_spec(6, extra_latex_after = "\\cmidrule(lr){2-11}") %>%
  row_spec(9, extra_latex_after = "\\cmidrule(lr){2-11}")

# Make the table font smaller using kable_styling()
latex_table <- latex_table %>%
  kable_styling(latex_options = c("scale_down"), font_size = 8)

# Print the LaTeX table code
cat(latex_table)


###################################################################################################################


res_df %>%
  summarize(kappa_IM = mean(kappa_IM), .by = n)




