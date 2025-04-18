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
registerDoMC(cores = 3)
#registerDoMC(cores = availableCores()-1)
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
      target_args <- list(W = W,
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
                          dx = 1e-4,
                          max_iter = 2000,
                          seq = TRUE,
                          nlambda_max = 10,
                          verbose = FALSE,
                          browse = FALSE)
      target_args_relax <- target_args; target_args_relax$grad_proj <- TRUE; target_args_relax$method <- "relax analytic"
      target_args_proj <- target_args; target_args_proj$grad_proj <- TRUE; target_args_proj$method <- "weak reg tmle"
      target_args_delta <- target_args; target_args_delta$grad_proj <- FALSE; target_args_delta$method <- "weak reg tmle"
      tau_relax <- do.call(target, target_args_relax)
      tau_proj <- do.call(target, target_args_proj)
      tau_delta <- do.call(target, target_args_delta)

      # point estimate and inference
      res_df <- map_dfr(seq(length(tau_proj)), function(.j) {
        .relax <- tau_relax[[.j]]
        .proj <- tau_proj[[.j]]
        .delta <- tau_delta[[.j]]

        psi_relax <- mean(.relax$pred)
        psi_proj <- mean(.proj$pred)
        psi_delta <- mean(.delta$pred)
        #Using the non parametric eic
        eic_relax <- .relax$eic
        eic_proj <- .proj$eic
        eic_delta <- .delta$eic

        se_relax_np <- sqrt(var(eic_relax, na.rm = TRUE) / .n)
        se_proj_np <- sqrt(var(eic_proj, na.rm = TRUE) / .n)
        se_delta_np<- sqrt(var(eic_delta, na.rm = TRUE) / .n)

        lower_relax_np <- psi_relax - 1.96 * se_relax_np
        upper_relax_np <- psi_relax + 1.96 * se_relax_np
        lower_proj_np <- psi_proj - 1.96 * se_proj_np
        upper_proj_np <- psi_proj + 1.96 * se_proj_np
        lower_delta_np <- psi_delta - 1.96 * se_delta_np
        upper_delta_np <- psi_delta + 1.96 * se_delta_np

        #Using the projected eic (weak)
        eic_relax <- .relax$eic_proj
        eic_proj <- .proj$eic_proj
        eic_delta <- .delta$eic_proj

        se_relax_proj <- sqrt(var(eic_relax, na.rm = TRUE) / .n)
        se_proj_proj <- sqrt(var(eic_proj, na.rm = TRUE) / .n)
        se_delta_proj <- sqrt(var(eic_delta, na.rm = TRUE) / .n)

        lower_relax_proj <- psi_relax - 1.96 * se_relax_proj
        upper_relax_proj <- psi_relax + 1.96 * se_relax_proj
        lower_proj_proj <- psi_proj - 1.96 * se_proj_proj
        upper_proj_proj <- psi_proj + 1.96 * se_proj_proj
        lower_delta_proj <- psi_delta - 1.96 * se_delta_proj
        upper_delta_proj <- psi_delta + 1.96 * se_delta_proj

        #Using the projected eic (cv)
        eic_relax <- .relax$eic_proj_cv
        eic_proj <- .proj$eic_proj_cv
        eic_delta <- .delta$eic_proj_cv

        se_relax_proj_cv <- sqrt(var(eic_relax, na.rm = TRUE) / .n)
        se_proj_proj_cv <- sqrt(var(eic_proj, na.rm = TRUE) / .n)
        se_delta_proj_cv <- sqrt(var(eic_delta, na.rm = TRUE) / .n)

        lower_relax_proj_cv <- psi_relax - 1.96 * se_relax_proj_cv
        upper_relax_proj_cv <- psi_relax + 1.96 * se_relax_proj_cv
        lower_proj_proj_cv <- psi_proj - 1.96 * se_proj_proj_cv
        upper_proj_proj_cv <- psi_proj + 1.96 * se_proj_proj_cv
        lower_delta_proj_cv <- psi_delta - 1.96 * se_delta_proj_cv
        upper_delta_proj_cv <- psi_delta + 1.96 * se_delta_proj_cv

        #Using the delta eic
        eic_relax <- .relax$eic_delta
        eic_proj <- .proj$eic_delta
        eic_delta <- .delta$eic_delta

        se_relax_delta <- sqrt(var(eic_relax, na.rm = TRUE) / .n)
        se_proj_delta <- sqrt(var(eic_proj, na.rm = TRUE) / .n)
        se_delta_delta <- sqrt(var(eic_delta, na.rm = TRUE) / .n)

        lower_relax_delta <- psi_relax - 1.96 * se_relax_delta
        upper_relax_delta <- psi_relax + 1.96 * se_relax_delta
        lower_proj_delta <- psi_proj - 1.96 * se_proj_delta
        upper_proj_delta <- psi_proj + 1.96 * se_proj_delta
        lower_delta_delta <- psi_delta - 1.96 * se_delta_delta
        upper_delta_delta <- psi_delta + 1.96 * se_delta_delta

        return(data.frame(n = .n,
                          B = .b,
                          j = .j,
                          relax_lambda = as.numeric(.relax$lambda[1]),
                          proj_lambda = as.numeric(.proj$lambda[1]),
                          delta_lambda = as.numeric(.delta$lambda[1]),
                          psi_relax = psi_relax,
                          psi_proj = psi_proj,
                          psi_delta = psi_delta,
                          ###non parametric eic based inference
                          se_relax_np = se_relax_np,
                          se_proj_np = se_proj_np,
                          se_delta_np = se_delta_np,
                          lower_relax_np = lower_relax_np,
                          upper_relax_np = upper_relax_np,
                          lower_proj_np = lower_proj_np,
                          upper_proj_np = upper_proj_np,
                          lower_delta_np = lower_delta_np,
                          upper_delta_np = upper_delta_np,
                          ###projected eic based inference (weak)
                          se_relax_proj = se_relax_proj,
                          se_proj_proj = se_proj_proj,
                          se_delta_proj = se_delta_proj,
                          lower_relax_proj = lower_relax_proj,
                          upper_relax_proj = upper_relax_proj,
                          lower_proj_proj = lower_proj_proj,
                          upper_proj_proj = upper_proj_proj,
                          lower_delta_proj = lower_delta_proj,
                          upper_delta_proj = upper_delta_proj,
                          ###projected eic based inference (cv)
                          se_relax_proj_cv = se_relax_proj_cv,
                          se_proj_proj_cv = se_proj_proj_cv,
                          se_delta_proj_cv = se_delta_proj_cv,
                          lower_relax_proj_cv = lower_relax_proj_cv,
                          upper_relax_proj_cv = upper_relax_proj_cv,
                          lower_proj_proj_cv = lower_proj_proj_cv,
                          upper_proj_proj_cv = upper_proj_proj_cv,
                          lower_delta_proj_cv = lower_delta_proj_cv,
                          upper_delta_proj_cv = upper_delta_proj_cv,
                          ###Delta eic based inference
                          se_relax_delta = se_relax_delta,
                          se_proj_delta = se_proj_delta,
                          se_delta_delta = se_delta_delta,
                          lower_relax_delta = lower_relax_delta,
                          upper_relax_delta = upper_relax_delta,
                          lower_proj_delta = lower_proj_delta,
                          upper_proj_delta = upper_proj_delta,
                          lower_delta_delta = lower_delta_delta,
                          upper_delta_delta = upper_delta_delta))
      })
      return(res_df)
    })
  })

  return(res_df)
}
