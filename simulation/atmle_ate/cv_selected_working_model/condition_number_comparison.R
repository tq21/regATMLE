.libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.2",
            .libPaths()))
`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
library(purrr)
library(torch)
library(origami)
library(hal9001)
library(glmnet)
library(furrr)
library(doMC)
library(devtools)
load_all()
#registerDoMC(cores = 3)
registerDoMC(cores = availableCores()-1)
source("dgp/sim_data.R")
set.seed(123)

run <- function(dgp_num) {
  B <- 200
  n_seq <- seq(500, 2000, 500)

  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n = " %+% .n %+% ", run: " %+% .b)
      data <- do.call(eval("sim_data_" %+% dgp_num), list(n = .n))
      W <- data[, grep("W", colnames(data))]
      A <- data$A
      Y <- data$Y
      folds <- make_folds(n = .n, V = 5)

      # estimate P(A=1|W)
      g1W <- learn_g(W = W,
                     A = A,
                     method = "glm",
                     folds = folds,
                     g_bounds = c(0.01, 0.99),
                     cross_fit_nuisance = TRUE)

      # estimate E(Y|W)
      theta <- learn_theta(W = W,
                           Y = Y,
                           delta = rep(1, .n),
                           method = "glm",
                           folds = folds,
                           family = "gaussian",
                           theta_bounds = NULL,
                           cross_fit_nuisance = TRUE)

      # estimate CATE
      tau_A <- learn_tau_A(W = W,
                           A = A,
                           Y = Y,
                           theta = theta,
                           g1W = g1W,
                           delta = rep(1, .n),
                           v_folds = 5,
                           weights = rep(1, .n),
                           enumerate_basis_args = list(max_degree = 2,
                                                       smoothness_orders = 1),
                           browse = FALSE)
      tau_delta <- target(W = W,
                          A = A,
                          Y = Y,
                          theta = theta,
                          g1W = g1W,
                          delta = rep(1, .n),
                          pseudo_outcome = tau_A$pseudo_outcome,
                          pseudo_weights = tau_A$pesudo_weights,
                          X = tau_A$X,
                          basis_list = tau_A$basis_list,
                          X_hal = tau_A$X_hal,
                          fit = tau_A$fit,
                          grad_proj = FALSE,
                          method = "weak reg tmle",
                          dx = 1e-4,
                          max_iter = 2000,
                          seq = FALSE,
                          nlambda_max = 1,
                          verbose = FALSE,
                          browse = FALSE)

      # point estimate and inference
      res_df <- map_dfr(seq(length(tau_delta)), function(.j) {
        .delta <- tau_delta[[.j]]
        psi_delta <- mean(.delta$pred)
        eic_delta <- eic_ate_wm(x_basis = .delta$x_basis,
                                g1W = g1W,
                                A = A,
                                Y = Y,
                                theta = theta,
                                tau = .delta$pred,
                                eic_method = "diag")
        kappa_IM <- eic_delta$kappa
        eic_delta <- eic_delta$eic
        se_delta <- sqrt(var(eic_delta, na.rm = TRUE) / .n)
        lower_delta <- psi_delta - 1.96 * se_delta
        upper_delta <- psi_delta + 1.96 * se_delta

        return(data.frame(n = .n,
                          B = .b,
                          j = .j,
                          delta_lambda = .delta$lambda,
                          kappa_IM = kappa_IM,
                          psi_delta = psi_delta,
                          se_delta = se_delta,
                          lower_delta = lower_delta,
                          upper_delta = upper_delta))
      })
      return(res_df)
    })
  })

  return(res_df)
}

res_df <- run(2)
write.csv(res_df, file = "out/res_kappa_dgp_" %+% 2 %+% "_" %+% timestamp %+% ".csv", row.names = FALSE)
