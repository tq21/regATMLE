library(purrr)
library(torch)
library(origami)
library(hal9001)
library(glmnet)
library(furrr)
library(doMC)
library(devtools)
load_all()
registerDoMC(cores = availableCores()-1)

source("sim_data.R")
set.seed(123)
B <- 50
n <- 500

res_df <- map_dfr(seq(B), function(.b) {
  print(.b)
  data <- sim_data(n)
  W <- data[, grep("W", colnames(data))]
  A <- data$A
  Y <- data$Y
  folds <- make_folds(n = n, V = 5)

  # estimate P(A=1|W)
  g1W <- learn_g(W = W,
                 A = A,
                 method = "glm",
                 folds = folds,
                 g_bounds = c(0.001, 0.999),
                 cross_fit_nuisance = TRUE)

  # estimate E(Y|W)
  theta <- learn_theta(W = W,
                       Y = Y,
                       delta = rep(1, n),
                       method = "glm",
                       folds = folds,
                       family = "gaussian",
                       theta_bounds = NULL,
                       cross_fit_nuisance = TRUE)

  # estimate CATE
  relax_an_args_list <- list(W = W,
                             A = A,
                             Y = Y,
                             theta = theta,
                             g1W = g1W,
                             delta = rep(1, n),
                             method = "relax analytic",
                             v_folds = 5,
                             weights = rep(1, n),
                             enumerate_basis_args = list(max_degree = 2,
                                                         smoothness_orders = 1),
                             verbose = FALSE)
  relax_gd_args_list <- relax_an_args_list
  relax_gd_args_list$method <- "relax gd"
  tau_an <- do.call(learn_tau_A, relax_an_args_list)
  tau_gd <- do.call(learn_tau_A, relax_gd_args_list)

  # point estimate and inference
  psi_relax_an <- mean(tau_an$pred)
  psi_relax_gd <- mean(tau_gd$pred)
  eic_relax_an_args_list <- list(tau_A = tau_an,
                                 g1W = g1W,
                                 theta = theta,
                                 Y = Y,
                                 A = A,
                                 eic_method = "diag")
  eic_relax_gd_args_list <- eic_relax_an_args_list
  eic_relax_gd_args_list$tau_A <- tau_gd
  eic_relax_an <- do.call(get_atmle_eic_psi, eic_relax_an_args_list)
  eic_relax_gd <- do.call(get_atmle_eic_psi, eic_relax_gd_args_list)
  PnEIC_relax_an <- mean(eic_relax_an)
  PnEIC_relax_gd <- mean(eic_relax_gd)
  se_relax_an <- sqrt(var(eic_relax_an, na.rm = TRUE) / n)
  se_relax_gd <- sqrt(var(eic_relax_gd, na.rm = TRUE) / n)
  lower_relax_an <- psi_relax_an - 1.96 * se_relax_an
  upper_relax_an <- psi_relax_an + 1.96 * se_relax_an
  lower_relax_gd <- psi_relax_gd - 1.96 * se_relax_gd
  upper_relax_gd <- psi_relax_gd + 1.96 * se_relax_gd

  return(data.frame(B = .b,
                    psi_relax_an = psi_relax_an,
                    psi_relax_gd = psi_relax_gd,
                    PnEIC_relax_an = PnEIC_relax_an,
                    PnEIC_relax_gd = PnEIC_relax_gd,
                    se_relax_an = se_relax_an,
                    se_relax_gd = se_relax_gd,
                    lower_relax_an = lower_relax_an,
                    upper_relax_an = upper_relax_an,
                    lower_relax_gd = lower_relax_gd,
                    upper_relax_gd = upper_relax_gd))
})
