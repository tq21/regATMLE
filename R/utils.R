.bound <- function(x, bounds) {
  return(pmin(pmax(x, bounds[1]), bounds[2]))
}

#' @title Make counterfactual design matrix
#'
#' @keywords internal
#'
#' @importFrom hal9001 make_design_matrix
#' @importFrom purrr is_empty
#'
#' @param basis_list A \code{list} of HAL bases.
#' @param X_counterfactual A \code{matrix} containing baseline covariates
#' evaluated at the desired counterfactual.
#' @param X_unpenalized A \code{matrix} containing baseline covariates or
#' bases that are not penalized by HAL. Default is \code{NULL}.
make_counter_design_matrix <- function(basis_list, X_counterfactual, X_unpenalized) {
  if (!is.null(X_unpenalized)) {
    if (is_empty(basis_list)) {
      # intercept + main terms
      return(cbind(1, X_unpenalized))
    } else {
      # intercept + HAL bases + main terms
      return(cbind(1, cbind(as.matrix(make_design_matrix(X_counterfactual, basis_list)), X_unpenalized)))
    }
  } else {
    if (is_empty(basis_list)) {
      # intercept
      return(matrix(rep(1, nrow(X_counterfactual))))
    } else {
      # intercept + HAL bases
      return(cbind(1, as.matrix(make_design_matrix(X_counterfactual, basis_list))))
    }
  }
}

svd_pseudo_inv <- function(mat, tol = 1e-6) {
  svd_res <- svd(mat)
  D_inv <- ifelse(svd_res$d > tol, 1 / svd_res$d, 0)
  return(svd_res$v %*% diag(D_inv) %*% t(svd_res$u))
}
