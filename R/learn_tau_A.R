#' @title Learn working model for the conditional effect of treatment on the
#' outcome
#'
#' @description Function to learn the conditional effect of treatment on the
#' outcome, \eqn{T(W)=\mathbb{E}(Y\mid W,A=1)-\mathbb{E}(Y\mid W,A=0)}.
#'
#' @keywords working model
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table
#'
#' @param W A matrix of baseline covariates.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param Y A vector of outcomes.
#' @param g A vector of estimated treatment probabilities,
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}.
#' @param theta_tilde A vector of estimated conditional mean of outcome given
#' baseline covariates, \eqn{\theta(W)=\mathbb{E}(Y\mid W)}.
#' @param method Working model type. Either \code{"glmnet"} for lasso-based
#' working model or \code{"HAL"} for highly adaptive lasso-based working model.
#' Default is \code{"glmnet"}.
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
#'
#' @returns A \code{list} containing the following elements:
#' \item{pred}{A numeric vector of the estimated conditional effects;}
#' \item{x_basis}{A numeric matrix of the working model bases;}
#' \item{coefs}{A numeric vector of the working model coefficients.}
learn_tau_A <- function(W,
                        A,
                        Y,
                        theta,
                        g1W,
                        delta,
                        v_folds,
                        weights,
                        enumerate_basis_args,
                        browse = FALSE) {

  if (browse) browser()

  # R-transformations
  pseudo_outcome <- ifelse(abs(A - g1W) < 1e-10, 0, (Y - theta) / (A - g1W))
  pseudo_weights <- (A - g1W)^2 * weights

  # check arguments
  enumerate_basis_default_args <- list(
    max_degree = 2,
    smoothness_orders = 1,
    num_knots = 20
  )
  enumerate_basis_args <- modifyList(
    enumerate_basis_default_args,
    enumerate_basis_args
  )

  X <- data.frame(W)

  # make design matrix
  basis_list <- enumerate_basis(x = as.matrix(X[delta == 1, ]),
                                max_degree = enumerate_basis_args$max_degree,
                                smoothness_orders = enumerate_basis_args$smoothness_orders,
                                num_knots = enumerate_basis_args$num_knots)
  X_hal <- make_design_matrix(X = as.matrix(X[delta == 1, ]),
                              blist = basis_list)

  # fit penalized HAL
  fit <- cv.glmnet(x = as.matrix(X_hal),
                   y = pseudo_outcome[delta == 1],
                   weights = pseudo_weights[delta == 1],
                   family = "gaussian",
                   alpha = 1,
                   nfolds = v_folds,
                   parallel = TRUE)

  return(list(pseudo_outcome = pseudo_outcome,
              pseudo_weights = pseudo_weights,
              X = X,
              basis_list = basis_list,
              X_hal = X_hal,
              fit = fit))
}
