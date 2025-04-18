glm_torch <- function(X,
                      Y,
                      family,
                      beta_init,
                      lr = 1e-4,
                      max_iter = 50000,
                      patience = 10,
                      tolerance = 1e-6,
                      verbose = TRUE) {

  # convert to torch tensors
  X <- torch_tensor(as.matrix(X))
  Y <- torch_tensor(Y)
  intercept <- torch_ones(X$size()[1], 1)
  phi_X <- torch_cat(list(intercept, X), dim = 2L)

  # determine loss function
  if (family == "gaussian") {
    loss_fn <- function(pred, Y) {
      return(pred$sub(Y)$pow(2)$mean())
    }
  } else if (family == "binomial") {
    loss_fn <- nn_bce_loss()
  }

  # define optimizer and LR scheduler
  beta <- torch_tensor(beta_init, requires_grad = TRUE)
  beta_optim <- optim_adam(list(beta), lr = lr)
  beta_scheduler <- lr_reduce_on_plateau(beta_optim, mode = 'min',
                                         patience = patience, factor = 0.5,
                                         verbose = verbose)

  # gradient descent
  train_losses <- numeric(max_iter)
  no_improve_counter <- 0
  best_loss <- Inf
  for (j in 1:max_iter) {
    beta_optim$zero_grad()
    pred <- phi_X$matmul(beta)
    train_loss <- loss_fn(pred, Y)
    train_loss$backward()
    beta_optim$step()
    train_losses[j] <- train_loss$item()

    if (train_loss$item() < best_loss - tolerance) {
      best_loss <- train_loss$item()
      no_improve_counter <- 0
    } else {
      no_improve_counter <- no_improve_counter + 1
    }

    if (no_improve_counter >= patience) {
      if (verbose) {
        cat("Early stopping at iteration", j, "with train loss:", train_loss$item(), "\n")
      }
      train_losses <- train_losses[1:j]
      break
    }

    if (verbose && j %% 100 == 0) {
      cat("iteration", j, "train loss:", train_loss$item(), "\n")
    }

    beta_scheduler$step(train_loss)
  }

  return(list(beta = as.numeric(beta)))
}
