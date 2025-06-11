plateau_selector_ci_based <- function(res_df,
                                      sort_by_l1_norm = FALSE) {
  # TODO: put a cap on the max se allowed
  # whether to sort by l1 norm or lambda
  if (sort_by_l1_norm) {
    res_df <- arrange(res_df, l1_norm)
    lambda <- "l1_norm"
  } else {
    lambda <- "lambda"
  }

  # estimate trend of psi
  if (nrow(res_df) == 1) {
    return(res_df)
  }
  fit_trend <- glm.fit(x = as.matrix(-res_df[[lambda]]),
                       y = res_df$psi,
                       family = gaussian())
  if (as.numeric(coef(fit_trend)) >= 0) {
    # maximize lower bound
    opt_idx <- which.max(res_df$lower)
  } else {
    # minimize upper bound
    opt_idx <- which.min(res_df$upper)
  }
  return(res_df[opt_idx,,drop=FALSE])
}
