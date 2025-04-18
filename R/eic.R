get_atmle_eic_psi <- function(tau_A,
                              g1W,
                              theta,
                              Y,
                              A,
                              eic_method = "diag") {
  IM <- t(tau_A$x_basis) %*% diag(g1W*(1-g1W)) %*% tau_A$x_basis / nrow(tau_A$x_basis)
  if (dim(tau_A$x_basis)[2] == 1) {
    IM_inv <- solve(IM)
  } else {
    if (eic_method == "svd_pseudo_inv") {
      # SVD-based pseudo-inverse
      IM_inv <- svd_pseudo_inv(IM)
    } else if (eic_method == "diag") {
      IM_inv <- solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
    }
  }
  D_beta <- as.vector(tau_A$x_basis %*% IM_inv %*% colMeans(tau_A$x_basis) * (A-g1W)*(Y-theta-(A-g1W)*tau_A$pred))
  W_comp <- tau_A$pred - mean(tau_A$pred)

  return(list(eic = W_comp + D_beta,
              kappa = kappa(IM)))
}

eic_ate <- function(QW1,
                    QW0,
                    psi,
                    A,
                    g1W,
                    Y,
                    QWA) {
  W_comp <- QW1-QW0-psi
  Y_comp <- (A/g1W-(1-A)/(1-g1W))*(Y-QWA)
  return(W_comp + Y_comp)
}

eic_ate_wm <- function(x_basis,
                       g1W,
                       A,
                       Y,
                       theta,
                       tau,
                       eic_method = "diag") {

  IM <- t(x_basis) %*% diag(g1W * (1 - g1W)) %*% x_basis / nrow(x_basis)
  if (dim(x_basis)[2] == 1) {
    IM_inv <- solve(IM)
  } else {
    if (eic_method == "svd_pseudo_inv") {
      # SVD-based pseudo-inverse
      IM_inv <- svd_pseudo_inv(IM)
    } else if (eic_method == "diag") {
      IM_inv <- solve(IM + diag(1e-6, nrow(IM), ncol(IM)))
    }
  }
  D_beta <- as.vector(x_basis %*% IM_inv %*% colMeans(x_basis)*(A-g1W)*(Y-theta-(A-g1W)*tau))
  W_comp <- tau - mean(tau)

  return(list(eic = W_comp + D_beta,
              kappa = kappa(IM)))
}


eic_ate_wm_TMLE <- function(x_basis,
                            QWA,
                            QW1,
                            QW0,
                            phi_0W_cur,
                            phi_1W_cur,
                            Y,
                            eic_method = "diag") {

  IM <- t(2*(Y-QWA)*x_basis)%*%(2*(Y-QWA)*x_basis)/nrow(x_basis)
  if (dim(x_basis)[2] == 1) {
    IM_inv <- solve(IM)
  } else {
    if (eic_method == "svd_pseudo_inv") {
      # SVD-based pseudo-inverse
      IM_inv <- svd_pseudo_inv(IM)
    } else if (eic_method == "diag") {
      IM_inv <- solve(IM + diag(1e-6, nrow(IM), ncol(IM)))
    }
  }


  D_beta <- as.vector(  t(as.vector(colMeans(as.matrix(phi_1W_cur[, -1, drop = FALSE]-phi_0W_cur[, -1, drop = FALSE])))) %*% IM_inv %*% t(2*(Y-QWA)*x_basis) )
  W_comp <- QW1-QW0-mean(QW1-QW0)

  dir<-as.vector( t(as.vector(colMeans(as.matrix(phi_1W_cur[, -1, drop = FALSE]-phi_0W_cur[, -1, drop = FALSE])))) %*% IM_inv)
  return(list(eic = W_comp + D_beta,
              kappa = kappa(IM),
              dir=dir))
}
