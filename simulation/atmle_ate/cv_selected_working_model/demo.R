.libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.2",
            .libPaths()))
library(purrr)
library(torch)
library(origami)
library(hal9001)
library(glmnet)
library(furrr)
library(doMC)
library(ggplot2)
library(devtools)
load_all()
registerDoMC(cores = availableCores()-1)
`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
source("sim_data.R")

#set.seed(123)
set.seed(93284)
n <- 500
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
               #g_bounds = c(0.001, 0.999),
               g_bounds = c(0, 1),
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
tau_A <- learn_tau_A(W = W,
                     A = A,
                     Y = Y,
                     theta = theta,
                     g1W = g1W,
                     delta = rep(1, n),
                     v_folds = 5,
                     weights = rep(1, n),
                     enumerate_basis_args = list(max_degree = 2,
                                                 smoothness_orders = 1),
                     browse = FALSE)
target_args <- list(W = W,
                    A = A,
                    Y = Y,
                    theta = theta,
                    g1W = g1W,
                    delta = rep(1, n),
                    pseudo_outcome = tau_A$pseudo_outcome,
                    pseudo_weights = tau_A$pesudo_weights,
                    X = tau_A$X,
                    basis_list = tau_A$basis_list,
                    X_hal = tau_A$X_hal,
                    fit = tau_A$fit,
                    dx = 1e-4,
                    max_iter = 2000,
                    seq = TRUE,
                    nlambda_max = 100,
                    verbose = FALSE,
                    browse = FALSE)
target_args_proj <- target_args; target_args_proj$grad_proj <- TRUE; target_args_proj$method <- "weak reg tmle"
target_args_delta <- target_args; target_args_delta$grad_proj <- FALSE; target_args_delta$method <- "weak reg tmle"
tau_proj <- do.call(target, target_args_proj)
tau_delta <- do.call(target, target_args_delta)

# point estimate and inference
res_df <- map_dfr(seq(length(tau_proj)), function(.j) {
  .proj <- tau_proj[[.j]]
  .delta <- tau_delta[[.j]]

  psi_proj <- mean(.proj$pred)
  psi_delta <- mean(.delta$pred)
  eic_proj <- eic_ate(QW1 = theta+(1-g1W)*.proj$pred,
                      QW0 = theta-g1W*.proj$pred,
                      psi = psi_proj,
                      A = A,
                      g1W = g1W,
                      Y = Y,
                      QWA = theta+(A-g1W)*.proj$pred)
  eic_delta <- eic_ate_wm(x_basis = .delta$x_basis,
                          g1W = g1W,
                          A = A,
                          Y = Y,
                          theta = theta,
                          tau = .delta$pred,
                          eic_method = "diag")
  se_proj <- sqrt(var(eic_proj, na.rm = TRUE) / n)
  se_delta <- sqrt(var(eic_delta, na.rm = TRUE) / n)
  lower_proj <- psi_proj - 1.96 * se_proj
  upper_proj <- psi_proj + 1.96 * se_proj
  lower_delta <- psi_delta - 1.96 * se_delta
  upper_delta <- psi_delta + 1.96 * se_delta

  return(data.frame(j = .j,
                    proj_lambda = .proj$lambda,
                    delta_lambda = .delta$lambda,
                    psi_proj = psi_proj,
                    psi_delta = psi_delta,
                    se_proj = se_proj,
                    se_delta = se_delta,
                    lower_proj = lower_proj,
                    upper_proj = upper_proj,
                    lower_delta = lower_delta,
                    upper_delta = upper_delta))
})

ggplot(res_df, aes(x = j, y = psi_proj)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_proj, ymax = upper_proj), width = 0.2) +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 2.5)) +
  labs(title = "Projection", x = "j", y = "ATE")

ggplot(res_df, aes(x = j, y = psi_delta)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_delta, ymax = upper_delta), width = 0.2) +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 2.5)) +
  labs(title = "Delta-method", x = "j", y = "ATE")
