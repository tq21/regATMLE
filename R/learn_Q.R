#' @title Learn nuisance function: conditional mean of outcome given treatment
#' and baseline covariates
#'
#' @description Function to learn the conditional mean of outcome given
#' baseline covariates, \eqn{Q(A,W)=\mathbb{E}(Y\mid A,W)}.
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
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param Y A vector of outcomes.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' @param folds A `fold` object from `origami` package.
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}. Either \code{"gaussian"} or \code{"binomial"}.
#' @param Q_bounds A numeric vector of lower and upper bounds for the
#' conditional mean of outcome given treatment and baseline covariates.
#' The first element is the lower bound, and the second element is the upper
#' bound.
#' @param cross_fit_nuisance A logical indicating whether to use cross-fitting
#' for estimating the nuisance function.
#'
#' @returns A `list` containing two numeric vectors:
#' \itemize{
#'   \item `A0`: Estimated conditional mean of outcome given control-arm
#'   \item `A1`: Estimated conditional mean of outcome given treatment-arm
#' }
learn_Q <- function(W,
                    A,
                    Y,
                    method,
                    folds,
                    family,
                    Q_bounds,
                    cross_fit_nuisance) {
  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners(family)
  }

  A0 <- A1 <- numeric(length(Y))

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    if (family == "gaussian") {
      lrnr <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                           Lrnr_cv_selector$new(loss_squared_error))
      task <- sl3_Task$new(data = data.table(W, A = A, Y = Y),
                           covariates = c(colnames(W), "A"), outcome = "Y",
                           folds = folds,
                           outcome_type = "continuous")
      task_A0 <- sl3_Task$new(data = data.table(W, A = 0, Y = Y),
                              covariates = c(colnames(W), "A"), outcome = "Y",
                              folds = folds,
                              outcome_type = "continuous")
      task_A1 <- sl3_Task$new(data = data.table(W, A = 1, Y = Y),
                              covariates = c(colnames(W), "A"), outcome = "Y",
                              folds = folds,
                              outcome_type = "continuous")
    } else if (family == "binomial") {
      lrnr <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                           Lrnr_cv_selector$new(loss_loglik_binomial))
      task <- sl3_Task$new(data = data.table(W, A = A, Y = Y),
                           covariates = c(colnames(W), "A"), outcome = "Y",
                           folds = folds,
                           outcome_type = "binomial")
      task_A0 <- sl3_Task$new(data = data.table(W, A = 0, Y = Y),
                              covariates = c(colnames(W), "A"), outcome = "Y",
                              folds = folds,
                              outcome_type = "binomial")
      task_A1 <- sl3_Task$new(data = data.table(W, A = 1, Y = Y),
                              covariates = c(colnames(W), "A"), outcome = "Y",
                              folds = folds,
                              outcome_type = "binomial")
    }
    fit <- lrnr$train(task)
    A0 <- .bound(fit$predict(task_A0), Q_bounds)
    A1 <- .bound(fit$predict(task_A1), Q_bounds)
  } else if (method == "glm") {
    X <- data.frame(W, A = A)
    X_A0 <- data.frame(W, A = 0)
    X_A1 <- data.frame(W, A = 1)

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(Y[train_idx] ~ ., data = X[train_idx, ], family = family)
        A0[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newdata = X_A0[valid_idx, ], type = "response"
        )), Q_bounds)
        A1[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newdata = X_A1[valid_idx, ], type = "response"
        )), Q_bounds)
      })
    } else {
      # no cross fit
      fit <- glm(Y ~ ., data = X, family = family)
      A0 <- .bound(as.numeric(predict(fit, newdata = X_A0, type = "response")),
                   Q_bounds)
      A1 <- .bound(as.numeric(predict(fit, newdata = X_A1, type = "response")),
                   Q_bounds)
    }
  } else {
    stop("Invalid method. Must be one of 'glm', or 'sl3', or a list of sl3
         learners.")
  }

  return(list(A0 = A0,
              A1 = A1))
}
