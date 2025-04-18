target_TMLE <- function(W,
                   A,
                   Y,
                   g1W,
                   QAW,
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

  #if (browse) browser()

  if (seq) {
    # consider a sequence of nested working models
    cv_lambda <- fit$lambda.min
    lambda_seq <- fit$lambda
    lambda_seq <- lambda_seq[lambda_seq <= cv_lambda]
    lambda_seq <- lambda_seq[1:min(nlambda_max, length(lambda_seq))] #lambda_seq[c(1,min(nlambda_max, length(lambda_seq)))]

    if(length(lambda_seq)==1){
      lambda_seq<-c(lambda_seq,lambda_seq)
    }

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
    phi_AW_cur <- cbind(1, X_hal_selected_cur)

    X0_hal<-QAW$X0_hal
    X0_hal_selected_cur <- X0_hal[, non_zero_cur, drop = FALSE]
    phi_0W_cur <- cbind(1, X0_hal_selected_cur)

    X1_hal<-QAW$X1_hal
    X1_hal_selected_cur <- X1_hal[, non_zero_cur, drop = FALSE]
    phi_1W_cur <- cbind(1, X1_hal_selected_cur)

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
        while (cur_iter <= max_iter & abs(PnEIC) > sn) {
          #if (browse) browser()
           # compute canonical gradient
          QW1 <- as.numeric(phi_1W_cur%*%beta_cur)
          QW0 <- as.numeric(phi_0W_cur%*%beta_cur)
          QWA <- as.numeric(phi_AW_cur%*%beta_cur)


          if (grad_proj) {
            # use the projection of the canonical gradient for TMLE update
            eic <- eic_ate(QW1 = QW1,
                           QW0 = QW0,
                           psi = mean(QW1-QW0),
                           A = A,
                           g1W = g1W,
                           Y = Y,
                           QWA = QWA)
            PnEIC <- mean(eic)
            score_mat <- (Y-QWA)*phi_AW_cur[,-1]
            proj_fit <- glmnet(x = score_mat,
                               y = eic,
                               intercept = FALSE,
                               lambda = 1e-5,
                               alpha = 1)
            direction <- as.numeric(coef(proj_fit, s = 1e-5)[-1])
          } else {
            # use the delta-method based gradient for TMLE update
            eic_list <- eic_ate_wm_TMLE(x_basis = as.matrix(phi_AW_cur[,-1]),
                                   QWA,
                                   QW1,
                                   QW0,
                                   phi_0W_cur,
                                   phi_1W_cur,
                                   Y,
                                   eic_method = "diag")
            eic<-eic_list$eic
            PnEIC <- mean(eic)
            direction <- eic_list$dir
          }

          # update beta
          ##################################################
          #Why we have the sign(PnEIC) here???????
          ##################################################
          beta_cur[-1] <- beta_cur[-1] + dx*sign(PnEIC)*direction
          sn <- sqrt(var(eic, na.rm = TRUE))/(sqrt(length(A)) * log(length(A)))
          cur_iter <- cur_iter + 1
          #print(PnEIC)
        }
      } else if (method == "cv reg tmle") {
        # CV-selected regularized targeting
        cur_iter <- 1
        PnEIC <- Inf
        sn <- 0
        while (cur_iter <= max_iter & abs(PnEIC) > sn) {
          # compute canonical gradient
          QW1 <- as.numeric(phi_1W_cur%*%beta_cur)
          QW0 <- as.numeric(phi_0W_cur%*%beta_cur)
          QWA <- as.numeric(phi_AW_cur%*%beta_cur)


          if (grad_proj) {
            # use the projection of the canonical gradient for TMLE update
            eic <- eic_ate(QW1 = QW1,
                           QW0 = QW0,
                           psi = mean(QW1-QW0),
                           A = A,
                           g1W = g1W,
                           Y = Y,
                           QWA = QWA)
            PnEIC <- mean(eic)
            score_mat <- (Y-QWA)*phi_AW_cur[,-1]
            proj_fit <- cv.glmnet(x = score_mat,
                                  y = eic,
                                  intercept = FALSE,
                                  parallel = TRUE,
                                  alpha = 1)
            direction <- as.numeric(coef(proj_fit, s = "lambda.min")[-1])
          } else {
            # use the delta-method based gradient for TMLE update
            eic_list <- eic_ate_wm_TMLE(x_basis = as.matrix(phi_AW_cur[,-1]),
                                   QWA,
                                   QW1,
                                   QW0,
                                   phi_0W_cur,
                                   phi_1W_cur,
                                   Y,
                                   eic_method = "diag")
            eic<-eic_list$eic
            PnEIC <- mean(eic)
            direction <- eic_list$dir
          }

          # update beta
          beta_cur[-1] <- beta_cur[-1] + dx*sign(PnEIC)*direction
          sn <- sqrt(var(eic, na.rm = TRUE))/(sqrt(length(A)) * log(length(A)))
          cur_iter <- cur_iter + 1
          #print(PnEIC)
        }
      }
      beta_cur[is.na(beta_cur)] <- 0
    } else {
      beta_cur <- mean(pseudo_outcome[delta == 1])
    }

    #design matrix of A,W
    x_basis <- make_counter_design_matrix(basis_list = basis_list_cur,
                                          X_counterfactual = as.matrix(X),
                                          X_unpenalized = NULL)
    #Q(A,W)
    pred <- as.numeric(x_basis%*%matrix(beta_cur))

    #Report Back
    # compute canonical gradient
    QW1 <- as.numeric(phi_1W_cur%*%beta_cur)
    QW0 <- as.numeric(phi_0W_cur%*%beta_cur)
    QWA <- as.numeric(phi_AW_cur%*%beta_cur)


    # use the projection of the canonical gradient for TMLE update
    eic <- eic_ate(QW1 = QW1,
                     QW0 = QW0,
                     psi = mean(QW1-QW0),
                     A = A,
                     g1W = g1W,
                     Y = Y,
                     QWA = QWA)

    score_mat <- (Y-QWA)*phi_AW_cur[,-1]

    #CV
    proj_fit <- cv.glmnet(x = score_mat,
                          y = eic,
                          intercept = FALSE,
                          parallel = TRUE,
                          alpha = 1)
    eic_proj_cv <-predict(proj_fit, newx = score_mat, s = "lambda.min")

    #Weak Regularization
    proj_fit <- glmnet(x = score_mat,
                       y = eic,
                       intercept = FALSE,
                       lambda = 1e-5,
                       alpha = 1)

    eic_proj <- predict(proj_fit, newx = score_mat)


    # use the delta-method based gradient for TMLE update
    eic_delta <- eic_ate_wm_TMLE(x_basis = as.matrix(phi_AW_cur[,-1]),
                             QWA,
                             QW1,
                             QW0,
                             phi_0W_cur,
                             phi_1W_cur,
                             Y,
                             eic_method = "diag")$eic
    eic_delta_kappa <- eic_ate_wm_TMLE(x_basis = as.matrix(phi_AW_cur[,-1]),
                                 QWA,
                                 QW1,
                                 QW0,
                                 phi_0W_cur,
                                 phi_1W_cur,
                                 Y,
                                 eic_method = "diag")$kappa


    return(list(lambda = lambda_seq[.j],
                pred=mean(QW1-QW0),
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
