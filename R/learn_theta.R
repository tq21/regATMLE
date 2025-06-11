#' @title Learn nuisance function: conditional mean of outcome given baseline
#' covariates
#'
#' @description Function to learn the conditional mean of outcome given
#' baseline covariates, \eqn{\tilde{\theta}(W)=\mathbb{E}(Y\mid W)}.
#'
#' @keywords nuisance
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table
#' @importFrom sl3 Stack
#' @importFrom sl3 make_learner
#' @importFrom sl3 sl3_Task
#' @importFrom sl3 Pipeline
#' @importFrom sl3 Lrnr_cv
#' @importFrom sl3 Lrnr_cv_selector
#' @importFrom sl3 loss_loglik_binomial
#' @importFrom sl3 loss_squared_error
#' @importFrom purrr walk
#'
#' @param W A matrix of baseline covariates.
#' @param Y A vector of outcomes.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' \code{"glmnet"} for lasso, or a \code{list} of \code{sl3} learners for
#' super learner-based estimation.
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}. Either \code{"gaussian"} or \code{"binomial"}.
#' @param theta_bounds A numeric vector of lower and upper bounds for the
#' conditional mean of outcome given baseline covariates.
#' The first element is the lower bound, and the second element is the upper
#' bound.
#'
#' @returns A numeric vector of the estimated values.
learn_theta <- function(W,
                        Y,
                        delta,
                        method,
                        folds,
                        family,
                        theta_bounds,
                        cross_fit_nuisance,
                        enumerate_basis_args) {
  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners(family)
  }

  pred <- numeric(length(Y))

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    if (family == "gaussian") {
      lrnr_theta <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                                 Lrnr_cv_selector$new(loss_squared_error))
      task_theta <- sl3_Task$new(data = data.table(W, Y = Y),
                                 covariates = colnames(W), outcome = "Y",
                                 folds = folds,
                                 outcome_type = "continuous")
    } else if (family == "binomial") {
      lrnr_theta <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                                 Lrnr_cv_selector$new(loss_loglik_binomial))
      task_theta <- sl3_Task$new(data = data.table(W, Y = Y),
                                 covariates = colnames(W), outcome = "Y",
                                 folds = folds,
                                 outcome_type = "binomial")
    }

    fit_theta <- lrnr_theta$train(task_theta)
    pred <- .bound(fit_theta$predict(task_theta), theta_bounds)
  } else if (method == "glm") {
    X <- data.frame(W)

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(Y[train_idx] ~ ., data = X[train_idx, ], family = family)
        pred[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newdata = X[valid_idx, ], type = "response"
        )), theta_bounds)
      })
    } else {
      # no cross fit
      fit <- glm(Y ~ ., data = X, family = family)
      pred <- .bound(
        as.numeric(predict(fit, newdata = X, type = "response")),
        theta_bounds
      )
    }
  } else if (method == "glmnet") {
    X <- as.matrix(W)
    fit <- cv.glmnet(x = X, y = Y,
                     keep = TRUE, alpha = 1, foldid = folds2foldvec(folds),
                     family = family)
    if (cross_fit_nuisance) {
      lambda_min <- fit$lambda[which.min(fit$cvm[!is.na(colSums(fit$fit.preval))])]
      pred <- as.numeric(fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == lambda_min])
      pred <- .bound(pred, theta_bounds)
    } else {
      pred <- .bound(as.numeric(predict(
        fit,
        newx = X, s = "lambda.min", type = "response"
      )), theta_bounds)
    }
  } else if (method == "HAL") {
    basis_list <- enumerate_basis(x = as.matrix(W),
                                  max_degree = enumerate_basis_args$max_degree,
                                  smoothness_orders = enumerate_basis_args$smoothness_orders,
                                  num_knots = enumerate_basis_args$num_knots)
    hal_design <- make_design_matrix(X = as.matrix(W), blist = basis_list)
    fit <- cv.glmnet(x = hal_design, y = Y,
                     keep = TRUE, alpha = 1, foldid = folds2foldvec(folds),
                     family = family, parallel = TRUE)
    if (cross_fit_nuisance) {
      lambda_min <- fit$lambda[which.min(fit$cvm[!is.na(colSums(fit$fit.preval))])]
      pred <- as.numeric(fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == lambda_min])
      pred <- .bound(pred, theta_bounds)
    } else {
      pred <- .bound(as.numeric(predict(
        fit,
        newx = hal_design, s = "lambda.min", type = "response"
      )), theta_bounds)
    }
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(pred)
}
