library(hal9001)
library(glmnet)
tlearner <- function(W,
                     A,
                     Y,
                     enumerate_basis_args = list(max_degree = 2,
                                                 smoothness_orders = 1,
                                                 num_knots = 10)) {

  # estimate E(Y|W,A=0)
  fit_Q0W <- fit_hal(X = as.matrix(W[A == 0,,drop=FALSE]),
                     Y = Y[A == 0],
                     max_degree = enumerate_basis_args$max_degree,
                     smoothness_orders = enumerate_basis_args$smoothness_orders,
                     num_knots = enumerate_basis_args$num_knots,
                     family = "gaussian")
  Q0W <- predict(fit_Q0W, new_data = as.matrix(W))

  # estimate CATE
  fit_cate <- fit_hal(W[A == 0,,drop=FALSE],
                      Y = Y[A == 0],
                      offset = Q0W[A == 0],
                      max_degree = enumerate_basis_args$max_degree,
                      smoothness_orders = enumerate_basis_args$smoothness_orders,
                      num_knots = enumerate_basis_args$num_knots,
                      family = "gaussian")
  tau <- predict(fit_cate, new_data = as.matrix(W))

  # compute E(Y|W,A=1)
  Q1W <- Q0W + tau

  return(list(Q0W = Q0W,
              Q1W = Q1W,
              tau = tau))
}
