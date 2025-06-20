target <- function(W,
                   A,
                   Y,
                   theta,
                   g1W,
                   delta,
                   method,
                   pseudo_outcome,
                   pseudo_weights,
                   X,
                   basis_list,
                   X_hal,
                   fit,
                   eic_method = "svd_pseudo_inv",
                   seq = FALSE,
                   nlambda_max = 10,
                   verbose = FALSE,
                   browse = FALSE) {

  if (browse) browser()

  if (seq) {
    # candidate lambda values
    cv_lambda <- fit$lambda.min
    lambda_seq <- fit$lambda
    lambda_seq <- lambda_seq[lambda_seq <= cv_lambda]
    lambda_seq <- lambda_seq[1:min(nlambda_max, length(lambda_seq))]

    # extract a sequence of nested working models
    non_zero_cur <- which(as.numeric(coef(fit, s = cv_lambda))[-1] != 0)
    non_zero_all <- list()
    lambda_wm_seq <- c()
    non_zero_all[[1]] <- non_zero_cur
    lambda_wm_seq[1] <- lambda_seq[1]
    for (j in 2:length(lambda_seq)) {
      non_zero_next <- which(as.numeric(coef(fit, s = lambda_seq[j]))[-1] != 0)
      #non_zero_all <- c(non_zero_all, list(non_zero_next))
      #lambda_wm_seq <- c(lambda_wm_seq, lambda_seq[j])
      if (!all(non_zero_next %in% non_zero_cur)) {
        non_zero_cur <- sort(union(non_zero_cur, non_zero_next)) # ensure working models are nested
        non_zero_all <- c(non_zero_all, list(non_zero_cur))
        lambda_wm_seq <- c(lambda_wm_seq, lambda_seq[j])
      }
    }
  } else {
    # CV selected working model
    lambda_seq <- lambda_wm_seq <- fit$lambda.min
    non_zero_all <- vector("list", 1)
    non_zero_all[[1]] <- which(as.numeric(coef(fit, s = lambda_seq))[-1] != 0)
  }

  # start with the cv selected as initial fit
  intercept <- as.numeric(coef(fit, s = fit$lambda.min))[1]
  beta <- as.numeric(coef(fit, s = fit$lambda.min))[-1]

  # extract a list of working models, up to nlambda_max beyond CV selected one
  # perform targeting in each working model
  res_list <- map(seq_along(non_zero_all), function(.j) {
    non_zero_cur <- non_zero_all[[.j]]
    lambda_cur <- lambda_wm_seq[.j]
    basis_list_cur <- basis_list[non_zero_cur]
    X_hal_selected_cur <- X_hal[, non_zero_cur, drop = FALSE]
    phi_W_cur <- cbind(1, X_hal_selected_cur)
    beta_cur <- c(intercept, beta[non_zero_cur])

    if (length(non_zero_cur) > 0) {
      if (method == "relaxed") {
        # relaxed fit
        fit_relax <- glm(pseudo_outcome[delta == 1] ~ .,
                         family = "gaussian",
                         data = data.frame(as.matrix(X_hal_selected_cur)),
                         weights = pseudo_weights[delta == 1])
        beta_cur <- as.numeric(coef(fit_relax))
      } else if (method == "oneshot") {
        x_basis <- as.matrix(phi_W_cur)
        n <- nrow(x_basis)
        IM <- t(x_basis) %*% diag(g1W*(1-g1W)) %*% x_basis / n
        IM_inv <- tryCatch({
          solve(IM)
        }, error = function(e) {
          if (eic_method == "svd_pseudo_inv") {
            svd_pseudo_inv(IM)
          } else if (eic_method == "diag") {
            solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
          } else {
            stop("Unknown eic_method specified.")
          }
        })
        clever_cov <- as.vector(IM_inv %*% colMeans(x_basis))
        H <- (A-g1W)*as.vector(x_basis %*% clever_cov)
        tau <- as.numeric(x_basis %*% beta_cur)
        R <- Y-theta-(A-g1W)*tau
        epsilon <- sum(H*R)/sum(H*H)
        beta_cur <- beta_cur+epsilon*clever_cov
      }
    } else {
      beta_cur <- mean(pseudo_outcome[delta == 1])
    }
    beta_cur[is.na(beta_cur)] <- 0
    cate_pred <- as.numeric(phi_W_cur%*%beta_cur)

    return(list(lambda = lambda_cur,
                pred = cate_pred,
                x_basis = as.matrix(phi_W_cur),
                coefs = beta_cur,
                non_zero = non_zero_cur))
  })

  return(res_list)
}
