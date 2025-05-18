`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
.libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.4",
            .libPaths()))
library(purrr)
library(origami)
library(hal9001)
library(glmnet)
library(furrr)
library(doMC)
library(data.table)
library(sl3)
library(devtools)
load_all()
options(sl3.verbose = FALSE)
registerDoMC(cores = availableCores()-1)
set.seed(123)
B <- 500
n_seq <- seq(500, 2000, 500)

run <- function(sim_data, gamma = 0.5) {
  # make sl3 learners
  learner_list <- list(
    Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0),
    Lrnr_xgboost$new(max_depth = 5, nrounds = 20, verbose = 0),
    Lrnr_ranger$new(),
    Lrnr_earth$new(degree = 2),
    Lrnr_gam$new()
  )

  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n = " %+% .n %+% ", run: " %+% .b)
      data <- sim_data(.n, gamma = gamma)
      W <- data[, grep("W", colnames(data))]
      A <- data$A
      Y <- data$Y
      folds <- make_folds(n = .n, V = 5)

      # estimate P(A=1|W)
      g1W <- learn_g(W = W,
                     A = A,
                     method = learner_list,
                     folds = folds,
                     g_bounds = c(0.01, 0.99),
                     cross_fit_nuisance = FALSE)

      # estimate E(Y|W)
      theta <- learn_theta(W = W,
                           Y = Y,
                           delta = rep(1, .n),
                           method = learner_list,
                           folds = folds,
                           family = "gaussian",
                           theta_bounds = NULL,
                           cross_fit_nuisance = FALSE)

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
      target_args <- list(W = W,
                          A = A,
                          Y = Y,
                          theta = theta,
                          g1W = g1W,
                          delta = rep(1, .n),
                          pseudo_outcome = tau_A$pseudo_outcome,
                          pseudo_weights = tau_A$pseudo_weights,
                          X = tau_A$X,
                          basis_list = tau_A$basis_list,
                          X_hal = tau_A$X_hal,
                          fit = tau_A$fit,
                          dx = 1e-5,
                          max_iter = 5000,
                          verbose = FALSE,
                          browse = FALSE)

      # relaxed-HAL
      target_args_relax <- target_args
      target_args_relax$method <- "relaxed"
      tau_relax <- do.call(target, target_args_relax)

      # arguments for HAL-TMLE
      target_args_tmle <- target_args
      target_args_tmle$method <- "oneshot"
      tau_tmle <- do.call(target, target_args_tmle)

      # point estimate and inference
      res_df <- map_dfr(seq(length(tau_relax)), function(.j) {
        .relax <- tau_relax[[.j]]
        .tmle <- tau_tmle[[.j]]
        psi_relax <- mean(.relax$pred)
        psi_tmle <- mean(.tmle$pred)
        eic_relax <- eic_ate_wm(x_basis = .relax$x_basis,
                                g1W = g1W,
                                A = A,
                                Y = Y,
                                theta = theta,
                                tau = .relax$pred,
                                eic_method = "svd_pseudo_inv")
        eic_tmle <- eic_ate_wm(x_basis = .tmle$x_basis,
                               g1W = g1W,
                               A = A,
                               Y = Y,
                               theta = theta,
                               tau = .tmle$pred,
                               eic_method = "svd_pseudo_inv")
        se_relax <- sqrt(var(eic_relax$eic, na.rm = TRUE) / .n)
        se_tmle <- sqrt(var(eic_tmle$eic, na.rm = TRUE) / .n)
        lower_relax <- psi_relax - 1.96 * se_relax
        upper_relax <- psi_relax + 1.96 * se_relax
        lower_tmle <- psi_tmle - 1.96 * se_tmle
        upper_tmle <- psi_tmle + 1.96 * se_tmle

        return(data.frame(n = .n,
                          B = .b,
                          j = .j,
                          psi_relax = psi_relax,
                          psi_tmle = psi_tmle,
                          lower_relax = lower_relax,
                          lower_tmle = lower_tmle,
                          upper_relax = upper_relax,
                          upper_tmle = upper_tmle))
      })
      return(res_df)
    })
  })

  return(res_df)
}
