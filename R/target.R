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
                   dx,
                   max_iter,
                   grad_proj = FALSE,
                   seq = FALSE,
                   nlambda_max = 10,
                   verbose = FALSE,
                   browse = FALSE) {

  if (browse) browser()

  if (seq) {
    # consider a sequence of nested working models
    cv_lambda <- fit$lambda.min
    lambda_seq <- fit$lambda
    lambda_seq <- lambda_seq[lambda_seq <= cv_lambda]
    lambda_seq <- lambda_seq[c(1,min(nlambda_max, length(lambda_seq)))] #lambda_seq[1:min(nlambda_max, length(lambda_seq))]

    # extract a sequence of nested working models
    non_zero_cur <- which(as.numeric(coef(fit, s = cv_lambda))[-1] != 0)
    non_zero_all <- list()
    non_zero_all[[1]] <- non_zero_cur
    for (j in 2:length(lambda_seq)) {
      non_zero_next <- which(as.numeric(coef(fit, s = lambda_seq[j]))[-1] != 0)
      if (!all(non_zero_next %in% non_zero_cur)) {
        non_zero_cur <- sort(union(non_zero_cur, non_zero_next)) # ensure nested working models
        non_zero_all <- c(non_zero_all, list(non_zero_cur))
      }
    }
  } else {
    # CV selected working model
    lambda_seq <- fit$lambda.min
    non_zero_all <- vector("list", 1)
    non_zero_all[[1]] <- which(as.numeric(coef(fit, s = lambda_seq))[-1] != 0)
  }
  intercept <- as.numeric(coef(fit, s = fit$lambda.min))[1]
  beta <- as.numeric(coef(fit, s = fit$lambda.min))[-1]

  # extract a list of working models, up to nlambda_max beyond CV selected one
  # perform targeting in each working model
  res_list <- map(seq_along(non_zero_all), function(.j) {
    non_zero_cur <- non_zero_all[[.j]]
    basis_list_cur <- basis_list[non_zero_cur]
    X_hal_selected_cur <- X_hal[, non_zero_cur, drop = FALSE]
    phi_W_cur <- cbind(1, X_hal_selected_cur)
    beta_cur <- c(intercept, beta[non_zero_cur])

    if (length(non_zero_cur) > 0) {
      if (method == "relax analytic") {
        # relaxed, using analytic formula
        fit_relax <- glm(pseudo_outcome[delta == 1] ~ .,
                         family = "gaussian",
                         data = data.frame(as.matrix(X_hal_selected_cur)),
                         weights = pseudo_weights[delta == 1])
        beta_cur <- as.numeric(coef(fit_relax))
      } else if (method == "relax gd") {
        # relaxed, using gradient descent
        fit_relax <- glm_torch(X = as.matrix(X_hal_selected_cur),
                               Y = pseudo_outcome[delta == 1],
                               family = "gaussian",
                               beta_init = beta_cur,
                               verbose = verbose)
        beta_cur <- fit$beta
      } else if (method == "weak reg tmle") {
        # weakly regularized targeting
        cur_iter <- 1
        PnEIC <- Inf
        sn <- 0

        # compute current PnEIC
        cate_pred <- as.numeric(phi_W_cur%*%beta_cur)
        eic <- eic_ate_wm(x_basis = as.matrix(phi_W_cur),
                          g1W = g1W,
                          A = A,
                          Y = Y,
                          theta = theta,
                          tau = cate_pred)$eic
        PnEIC_prev <- mean(eic)
        prev_pos <- PnEIC_prev >= 0

        while (cur_iter <= max_iter & abs(PnEIC_prev) > sn) {
          # compute canonical gradient
          cate_pred <- as.numeric(phi_W_cur%*%beta_cur)
          eic <- eic_ate_wm(x_basis = as.matrix(phi_W_cur),
                            g1W = g1W,
                            A = A,
                            Y = Y,
                            theta = theta,
                            tau = cate_pred)$eic

          # if signs flip, reduce step size
          cur_pos <- PnEIC_prev >= 0
          if (cur_pos != prev_pos) {
            dx <- dx/10
          }
          prev_pos <- cur_pos

          # if PnEIC diff to small, increase step size
          PnEIC_cur <- mean(eic)
          direction <- get_beta_h(x_basis = as.matrix(phi_W_cur), g1W = g1W)
          if (abs(PnEIC_cur - PnEIC_prev) < 1e-6) {
            dx <- dx*10
          }
          PnEIC_prev <- PnEIC_cur

          # update beta
          beta_cur <- beta_cur + dx*sign(PnEIC_cur)*direction
          sn <- 0.001*sqrt(var(eic, na.rm = TRUE))/(sqrt(length(A)) * log(length(A)))
          cur_iter <- cur_iter + 1
          print(PnEIC_cur)
        }
      } else if (method == "cv reg tmle") {
        # CV-selected regularized targeting
        cur_iter <- 1
        PnEIC <- Inf
        sn <- 0
        while (cur_iter <= max_iter & abs(PnEIC) > sn) {
          # compute canonical gradient
          cate_pred <- as.numeric(phi_W_cur%*%beta_cur)

          if (grad_proj) {
            # use the projection of the canonical gradient for TMLE update
            eic <- eic_ate(QW1 = theta+(1-g1W)*cate_pred,
                           QW0 = theta-g1W*cate_pred,
                           psi = mean(cate_pred),
                           A = A,
                           g1W = g1W,
                           Y = Y,
                           QWA = theta+(A-g1W)*cate_pred)
            PnEIC <- mean(eic)
            score_mat <- (Y-theta-(A-g1W)*cate_pred)*(A-g1W)*phi_W_cur
            proj_fit <- cv.glmnet(x = score_mat,
                                  y = eic,
                                  intercept = FALSE,
                                  parallel = TRUE,
                                  alpha = 1)
            direction <- as.numeric(coef(proj_fit, s = "lambda.min")[-1])
          } else {
            # use the delta-method based gradient for TMLE update
            eic <- eic_ate_wm(x_basis = as.matrix(phi_W_cur),
                              g1W = g1W,
                              A = A,
                              Y = Y,
                              theta = theta,
                              tau = cate_pred)$eic
            PnEIC <- mean(eic)
            direction <- get_beta_h(x_basis = as.matrix(phi_W_cur), g1W = g1W)
          }

          # update beta
          beta_cur <- beta_cur + dx*sign(PnEIC)*direction
          sn <- sqrt(var(eic, na.rm = TRUE))/(sqrt(length(A)) * log(length(A)))
          cur_iter <- cur_iter + 1
          #print(PnEIC)
        }
      }
      beta_cur[is.na(beta_cur)] <- 0
    } else {
      beta_cur <- mean(pseudo_outcome[delta == 1])
    }

    x_basis <- make_counter_design_matrix(basis_list = basis_list_cur,
                                          X_counterfactual = as.matrix(X),
                                          X_unpenalized = NULL)
    pred <- as.numeric(x_basis%*%matrix(beta_cur))



    # compute canonical gradient
    cate_pred <- as.numeric(phi_W_cur%*%beta_cur)
    # use the projection of the canonical gradient for TMLE update
    eic <- eic_ate(QW1 = theta+(1-g1W)*cate_pred,
                   QW0 = theta-g1W*cate_pred,
                   psi = mean(cate_pred),
                   A = A,
                   g1W = g1W,
                   Y = Y,
                   QWA = theta+(A-g1W)*cate_pred)


    score_mat <- (Y-theta-(A-g1W)*cate_pred)*(A-g1W)*phi_W_cur

    # use the delta-method based gradient for TMLE update
    eic_delta <- eic_ate_wm(x_basis = as.matrix(phi_W_cur),
                            g1W = g1W,
                            A = A,
                            Y = Y,
                            theta = theta,
                            tau = cate_pred)$eic

    eic_delta_kappa <- eic_ate_wm(x_basis = as.matrix(phi_W_cur),
                                  g1W = g1W,
                                  A = A,
                                  Y = Y,
                                  theta = theta,
                                  tau = cate_pred)$kappa


    if(dim(score_mat)[2]<2){
      eic_proj<-eic
      eic_proj_cv<-eic
    }else{
      #Weak Regularization
      proj_fit <- glmnet(x = score_mat,
                         y = eic,
                         intercept = FALSE,
                         lambda = 1e-5,
                         alpha = 1)

      eic_proj <- predict(proj_fit, newx = score_mat)

      #CV

      proj_fit <- cv.glmnet(x = score_mat,
                            y = eic,
                            intercept = FALSE,
                            parallel = TRUE,
                            alpha = 1)
      eic_proj_cv <-predict(proj_fit, newx = score_mat, s = "lambda.min")
    }



    return(list(lambda = lambda_seq[.j],
                pred = pred,
                eic = eic,
                eic_proj = eic_proj,
                eic_proj_cv = eic_proj_cv,
                eic_delta = eic_delta,
                eic_delta_kappa=eic_delta_kappa,
                x_basis = x_basis,
                coefs = beta_cur,
                non_zero = non_zero_cur))
  })

  return(res_list)
}
