# ATTENTION: gamma = 0.5, 1, or 2 to adjust degree of overlap
sim_data <- function(n,
                     gamma = 0.5,
                     counter_A = NULL) {
  # error
  UY <- rnorm(n, 0, sqrt(0.5))

  # baseline covariates
  W1 <- round(runif(n, -1, 1), 2)
  W2 <- round(runif(n, -1, 1), 2)
  W3 <- round(runif(n, -1, 1), 2)
  W4 <- round(runif(n, -1, 1), 2)

  # treatment
  if (is.null(counter_A)) {
    A <- rbinom(n, 1, plogis(gamma*(W1+W2+W3+W4+sin(4*W1)+sin(4*W2)+sin(4*W3)+sin(4*W4))))
  } else {
    A <- rep(counter_A, n)
  }

  # outcome
  Y <- W1+abs(W2)+W3+abs(W4)+A*(1+W1+abs(W2)+cos(4*W3)+W4)+UY

  data <- data.frame(W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     W4 = W4,
                     A = A,
                     Y = Y)

  return(data)
}

get_truth <- function() {
  data_A1 <- sim_data(1e7, counter_A = 1)
  data_A0 <- sim_data(1e7, counter_A = 0)
  return(mean(data_A1$Y- data_A0$Y))
}
