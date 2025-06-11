library(dplyr)
source("dgps/sim_data.R")
truth <- get_truth()

res_df <- read.csv("out/res_gamma_5_0611_161800.csv")
res_df %>%
  filter(j == 1) %>%
  summarize(abs_bias_relax = abs(mean(psi_relax - truth)),
            abs_bias_tmle = abs(mean(psi_tmle - truth)),
            se_relax = sd(psi_relax),
            se_tmle = sd(psi_tmle),
            bias_se_relax = abs_bias_relax/se_relax,
            bias_se_tmle = abs_bias_tmle/se_tmle,
            mse_relax = mean((psi_relax - truth)^2),
            mse_tmle = mean((psi_tmle - truth)^2),
            cover_relax = mean(truth >= lower_relax & truth <= upper_relax),
            cover_tmle = mean(truth >= lower_tmle & truth <= upper_tmle),
            oracle_cover_relax = mean(truth >= psi_relax - 1.96 * sd(psi_relax) & truth <= psi_relax + 1.96 * sd(psi_relax)),
            oracle_cover_tmle = mean(truth >= psi_tmle - 1.96 * sd(psi_tmle) & truth <= psi_tmle + 1.96 * sd(psi_tmle)),
            .by = "n")




# TODO: add lambda
# cv-selected estimate
res_df <- res_df %>%
  select(n, B, j, lambda, psi_tmle, lower_tmle, upper_tmle) %>%
  rename(psi = psi_tmle, lower = lower_tmle, upper = upper_tmle)

n_seq <- seq(500, 2000, 500)

cv_select_df <- map_dfr(n_seq, function(.n) {
  map_dfr(seq(B), function(.b) {
    tmp_df <- res_df %>% filter(n == .n, B == .b)
    return(tmp_df[1,,drop=FALSE])
  })
})
cv_select_df$select_method <- "cv"

# CI-based selector selected estimate
ci_based_select_df <- map_dfr(n_seq, function(.n) {
  map_dfr(seq(B), function(.b) {
    tmp_df <- res_df %>% filter(n == .n, B == .b)
    return(plateau_selector_ci_based(tmp_df))
  })
})
ci_based_select_df$select_method <- "ci_based"

select_df <- rbind(cv_select_df,
                   ci_based_select_df)
select_df %>%
  summarize(abs_bias = abs(mean(psi-truth)),
            se = sd(psi),
            mse = mean((psi-truth)^2),
            coverage = mean(lower <= truth & upper >= truth),
            .by = c("n", "select_method")) %>%
  arrange(n)



